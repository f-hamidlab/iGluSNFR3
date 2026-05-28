# iGluSNFR3

Detection of synaptic activity in widefield single-plane fluorescence imaging [iGluSNFR3](https://www.nature.com/articles/s41592-023-01863-6).

This repository currently provides two maintained MATLAB entry pipelines:

- `multi_cell_activity_detection_pipeline.m`: single-condition pipeline (supports spontaneous and evoked datasets via config selection)
- `multi_cell_activity_detection_pipeline_matching.m`: evoked multi-condition pipeline with post hoc cluster matching

Use the first pipeline when you want to analyze one condition at a time. Use the second pipeline when one folder contains multiple evoked conditions that should be matched after detection.

The code processes Bio-Formats-compatible microscopy files, detects events at pixel/ROI level, clusters synchronous activity, and exports figures and `processed_data.mat` outputs.

## Scope and Supported Data

- Single-plane, single-channel fluorescence recordings
- iGluSNFR3 probe datasets
- File formats readable by Bio-Formats (for example `.cxd`, `.tif`, `.nd2`)
- Typical acquisition rates around 100 Hz (set in config)

## Requirements

### 1) MATLAB

Use a recent MATLAB release with these toolboxes available:

- Image Processing Toolbox
- Signal Processing Toolbox

### 2) Included third-party code (already in this repository)

- Bio-Formats MATLAB utilities in `Scripts/bfmatlab/` (including `bioformats_package.jar`)
- MLspike helper code in `Scripts/MLspike/`

No separate clone is required if you use this repository as-is.

## Installation

1. Clone or download this repository.
2. Open MATLAB.
3. Set your MATLAB current folder to the repository root (the folder containing `multi_cell_activity_detection_pipeline.m`).
4. Open one of the pipeline entry scripts and edit paths/configuration before running (details below).

## Quick Start

### A) Single-condition pipeline (spontaneous or evoked)

Use `multi_cell_activity_detection_pipeline.m`.

1. In the script, set:
	 - `ops.filedir` to your raw data root
	 - `ops.savedir` to your output root
2. Choose the config based on dataset type:
	 - `ops = config_spontaneous(ops);` for spontaneous activity datasets
	 - `ops = config_evoked(ops);` for single-condition evoked datasets
3. In the selected config file, set input matching fields as needed:
	 - `ops.fileformat` (for example `.cxd` or `.tif`)
	 - `ops.filename_regex` (regex filter for accepted filenames)
4. Confirm acquisition frame rate (`ops.fs`) in the config file you selected.
5. Run the script in MATLAB.

Behavior:

- Recursively scans subfolders under `ops.filedir`
- Creates per-file output folders under `ops.savedir`
- Writes `processed_data.mat`, figure files, and optional masks/maps

### B) Evoked multi-condition pipeline with matching

Use `multi_cell_activity_detection_pipeline_matching.m`.

1. In the script, set:
	 - `ops.filedir`
	 - `ops.savedir`
2. Choose/configure one multi-spec config file:
	 - `Config/config_evoked_multi_spec.m` for real experiments
	 - `Config/config_evoked_multi_spec_test_data.m` for bundled test examples
3. In the chosen config file:
	 - set input matching fields:
		 - `ops.fileformat`
		 - `ops.filename_regex` (default expects names like `Cell1_1.tif`)
	 - set `N` to the number of images/conditions
	 - customize per-image stimulus parameters (`first_stim`, `n_stim`, `stim_freq`, `len_spike`) as needed
	 - verify `ops.fs`
4. Run `multi_cell_activity_detection_pipeline_matching.m`.

Behavior:

- Applies condition-specific parameters based on image number parsed from filename
- If image index exceeds available condition configs, selection loops cyclically (e.g. 5->1, 6->2 for 4 configs)
- Runs event detection for each file
- Runs cluster matching at the end of each leaf folder (currently via `matching_clusters_with_image4`)

## Minimal Example Run

The example below runs the spontaneous pipeline on one TIFF recording.

Example folder layout:

```text
iGluSNFR3/
├── multi_cell_activity_detection_pipeline.m
├── Config/
├── Scripts/
├── originals/
│   └── demo_spont/
│       └── Cell1_1.tif
└── outputs/
```

Edit `multi_cell_activity_detection_pipeline.m`:

```matlab
ops.filedir = './originals/demo_spont/';
ops.savedir = './outputs/demo_spont/';
ops = config_spontaneous(ops);
% Set ops.fileformat and ops.filename_regex in Config/config_spontaneous.m
```

Run in MATLAB:

```matlab
run('multi_cell_activity_detection_pipeline.m')
```

Expected output after a successful run:

```text
outputs/
└── demo_spont/
	└── Cell1_1/
		├── processed_data.mat
		├── ROI_pxMap/                  (created; may be empty if plotting is off)
		└── Fig*.png                    (if visualization is enabled)
```

Notes:

- If no figures are created, check whether `ops.visualize` is set to `false` in your config.
- If your file is not detected, verify `ops.fileformat` and `ops.filename_regex` in your selected config.

## Key Configuration Files

- `Config/config_spontaneous.m`: defaults for spontaneous analysis
- `Config/config_evoked.m`: defaults for single-condition evoked analysis
- `Config/config_evoked_multi_spec.m`: per-condition evoked settings for multi-spec experiments
- `Config/config_evoked_multi_spec_test_data.m`: test-data variant of multi-spec settings

Important: these configs define file matching (`ops.fileformat`, `ops.filename_regex`) and convert many time-based parameters into frames internally. Always verify `ops.fs` matches your acquisition settings before running.

## Advanced Option: Create Your Own Config File

You can create a custom config file for your lab or experiment instead of editing the default configs.

This is useful when you want to keep the shipped config files unchanged and put your own experiment settings in a separate file.

### 1) Create a new config in `Config/`

1. Copy a starting template:
   - `Config/config_spontaneous.m` for spontaneous-like data
   - `Config/config_evoked.m` for single-condition evoked data
2. Save as a new file, for example `Config/config_my_experiment.m`.
3. Rename the function inside the file so it matches the filename:

```matlab
function ops = config_my_experiment(ops)
```

### 2) Use your custom config in the single-condition pipeline

In `multi_cell_activity_detection_pipeline.m`, replace the config call with your custom function:

```matlab
ops = config_my_experiment(ops);
```

### 3) Optional: Use a custom multi-spec config for matching pipeline

For multi-condition evoked workflows, copy and customize:

- `Config/config_evoked_multi_spec.m`

Then call your custom multi-spec config in `multi_cell_activity_detection_pipeline_matching.m`:

```matlab
ops_multi = config_my_multi_spec(ops);
ops = ops_multi{1};
```

In this example, the custom multi-spec config returns one settings set per condition, so `ops_multi{1}` selects the first one.

### 4) Recommended validation checks

- Confirm `ops.fs` matches acquisition frame rate.
- Confirm `ops.experiment_type` is set appropriately.
- Run one short recording first and inspect `processed_data.mat` and figures.

## Input Data Notes

- If `ops.use_binary_mask = false`, the pipeline performs automatic threshold-based segmentation.
- If `ops.use_binary_mask = true`, mask lookup is configurable from the selected config file in `Config/`.

### Binary Mask Lookup Configuration

Set these fields in your selected config file (for example `Config/config_evoked.m`, `Config/config_spontaneous.m`, or multi-spec config defaults):

```matlab
ops.binary_mask_location = "auto";
ops.binary_mask_pattern_filedir = 'MAX_Cell*_binary_*.tif';
ops.binary_mask_pattern_savedir = 'Cell*_*_binary.tif';
ops.binary_mask_pattern_savedir_parent_recursive = 'Cell*_*_binary.tif';
ops.binary_mask_group_size = 4;
```

Supported `ops.binary_mask_location` values:

- `"filedir"`: search only in the same folder as the source image (`ops.filename`) using `ops.binary_mask_pattern_filedir`
- `"savedir"`: search only in the per-file output folder (`ops.savedir`) using `ops.binary_mask_pattern_savedir`
- `"savedir_parent_recursive"`: search recursively under `parent(ops.savedir)` using `ops.binary_mask_pattern_savedir_parent_recursive`
- `"auto"`: search all locations above in that priority order

Savedir-specific mask selection rule:

- For masks found via `"savedir"` or `"savedir_parent_recursive"`, the script expects filenames like:
	- `CellX_num_binary.tif` (for example `Cell1_3_binary.tif`)
- Let `n` be the recording number from the image filename `CellX_n.tif`.
- Let `g = ops.binary_mask_group_size`.
- The target mask number is computed as:
	- $floor((n-1)/g) * g + 1$
- Only masks whose `num` matches that target are considered valid candidates.

Note for multi-condition configs:

- In `config_evoked_multi_spec*.m`, `ops.binary_mask_group_size` is synced to top-level `N` by default.

Examples:

- Image `Cell1_1.tif` with `g=4` -> target mask number $floor((1-1)/4)*4+1 = 1$ -> selects `Cell1_1_binary.tif`
- Image `Cell1_4.tif` with `g=4` -> target mask number $floor((4-1)/4)*4+1 = 1$ -> selects `Cell1_1_binary.tif`
- Image `Cell1_5.tif` with `g=4` -> target mask number $floor((5-1)/4)*4+1 = 5$ -> selects `Cell1_5_binary.tif`

Behavior when multiple files match after filtering:

- The script keeps the first match found in search order and emits a warning.

## Output Structure

For each processed recording, output folders typically contain:

- `processed_data.mat`: detected events, ROI/spike statistics, analysis outputs
- Figure files (format from `ops.fig_format`)
- `ROI_pxMap/` (if enabled)
- Optional ImageJ-compatible TIFF masks (depending on save flags)

## Performance Tips

- Set `ops.visualize = false` to reduce runtime and memory usage.
- Keep `ops.close_fig = true` for large batch runs.
- Use `ops.redo_detection = false` to skip files that already have `processed_data.mat`.

## Troubleshooting

- Bio-Formats/Java errors:
	- ensure the repository root is your MATLAB current folder
	- ensure `Scripts/bfmatlab/` is added (the pipeline does this automatically)
- No files processed:
	- check `ops.fileformat`
	- check `ops.filename_regex`
	- verify your input folder tree under `ops.filedir`
- Poor detection quality:
	- verify frame rate (`ops.fs`)
	- tune thresholds in the relevant config file (`m_thres`, `SNR_thres`, peak parameters)

## Citation

If you use this analysis pipeline, please cite:

- 

## Additional Documentation

Project wiki: https://github.com/f-hamidlab/iGluSNFR3/wiki/

