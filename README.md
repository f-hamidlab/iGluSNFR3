# iGluSNFR3 analysis pipeline

MATLAB pipeline for detecting and analysing synaptic activity in single-plane
widefield microscopy recordings acquired with the iGluSNFR3 probe. The current
entry point in this repository is
`multi_cell_activity_detection_pipeline_evoked_matching.m`, which is configured
for evoked responses and can match event clusters across recordings.

The pipeline accepts Bio-Formats-compatible image files. The supplied demo uses
`.cxd` files and binary masks.

## 1. System requirements

### Operating system

- Supported operating system: Windows, Linux, macOS
- Operating system version used for testing: **[FILL IN: Windows 11]**

The repository has not yet recorded a verified cross-platform compatibility
matrix. The pipeline uses MATLAB graphics, Java-based Bio-Formats, and
parallel processing, so these components must be available on the selected
operating system.

### Software

- MATLAB: R2024b
- Image Processing Toolbox: required for `strel`, `imtophat`, `regionprops`,
  `bwlabel`, `imextendedmax`, and related image operations.
- Signal Processing Toolbox: required for functions such as `findpeaks` and
  signal-processing utilities.
- Statistics and Machine Learning Toolbox: required for `pdist2`, `linkage`,
  `cluster`, and related clustering functions.
- Parallel Computing Toolbox: required for the `parfor` pixel-analysis loop. (Optional but highly recommended as it significantly reduce computation time. If parallel computation is not available, change the`parfor` in the main script to `for`)
- Java runtime supported by the selected MATLAB release: required by
  Bio-Formats.
- Bio-Formats MATLAB toolbox: included in `Scripts/bfmatlab/`. The source
  documentation referenced by the code is [Bio-Formats MATLAB documentation](https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html).
  Bundled/tested Bio-Formats version: v7.0.0
- [MLspike](https://github.com/MLspike/spikes): included in `Scripts/MLspike/`. The pipeline calls
  `tps_mlspikes` while initialising spike-train parameters, even when
  `ops.ST.option` is set to `"findpeaks"`. Bundled/tested MLspike version:
  commit 32fb84e Apr 9, 2020

The required MATLAB functions are included in this repository unless listed as
part of a MATLAB toolbox above. No Python installation is required.

### Hardware

- No special hardware is required to run the analysis on existing files.
- Data acquisition requires a compatible widefield fluorescence microscope,
  an iGluSNFR3-labelled preparation, and recordings in a supported image
  format. These are not required for software-only analysis.
- Recommended CPU, RAM, and storage for the supplied recordings:
  **[FILL IN after benchmarking]**

### Tested configuration

Complete this section after a clean run:

- MATLAB release: **24.2.0.3070828 (R2024b) Update 7**
- MATLAB toolbox versions: **Version 24.2 (R2024b)**
- Bio-Formats version: **v7.0.0**
- MLspike version or commit: 32fb84e Apr 9, 2020
- Operating system and version: **[FILL IN: Windows 11]**
- CPU/RAM: **[FILL IN: Intel Xeon Silver 4214 CPU @ 2.20 GHz, 24 physical cores/48 logical CPUs, 251 GiB RAM]**

## 2. Installation guide

### Instructions

1. Obtain this repository and open the `iGluSNFR3` directory in MATLAB.
2. Confirm that the required MATLAB toolboxes listed above are installed:

   ```matlab
   ver
   ```

   You would see an output similar to the following: **[FILL IN: replace output for Windows]**

   ```matlab
   -------------------------------------------------------------------------------------------------------------
   MATLAB Version: 24.2.0.3070828 (R2024b) Update 7
   MATLAB License Number: 860095
   Operating System: Linux 6.8.0-139-generic #139-Ubuntu SMP PREEMPT_DYNAMIC Sat Aug  1 03:52:05 UTC 2026 x86_64
   Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
   -------------------------------------------------------------------------------------------------------------
   MATLAB                                                Version 24.2        (R2024b)
   Image Processing Toolbox                              Version 24.2        (R2024b)
   Parallel Computing Toolbox                            Version 24.2        (R2024b)
   Signal Processing Toolbox                             Version 24.2        (R2024b)
   Statistics and Machine Learning Toolbox               Version 24.2        (R2024b)
   ```

   If you don't have the required toolbox, you can follow the instructions [here](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html) to install them.

No package manager or compilation step is required for the bundled MATLAB
scripts. Installation time on a normal desktop computer: **[FILL IN, excluding
download time]**.

## 3. Demo

The repository includes three `.cxd` demo recordings in
`Test_data/originals/` and corresponding processed results in
`Test_data/outputs/`. The supplied recordings are `Cell1_1.cxd`, `Cell1_2.cxd`,
and `Cell1_3.cxd`. The binary mask used by the pipeline is
`MAX_Cell1_1_binary_dil3.tif` in the same directory.

### Instructions to run the demo

1. Start MATLAB and change directory to the repository root.
2. Open `multi_cell_activity_detection_pipeline_evoked_matching.m`.
3. In the **Defining parameters and file paths** section, set:

   ```matlab
   % ========== FILE PATHS ==========
   % path to data ; The script will loop through all subfolders.
   scriptDir = fileparts(mfilename('fullpath')); % no need to change
   ops.filedir = fullfile(scriptDir, 'Test_data', 'originals');
   ops.fileformat = '.cxd';

   % path to saving directory
   ops.savedir = fullfile(scriptDir, 'Test_data', 'outputs');

   % path to other required functions
   addpath(fullfile(scriptDir, 'Scripts'))
   addpath(fullfile(scriptDir, 'Scripts/bfmatlab/'))
   addpath(fullfile(scriptDir, 'Scripts/MLspike/brick/'))
   addpath(fullfile(scriptDir, 'Scripts/MLspike/spikes/'))
   ```
4. Run the script. It recursively processes the `.cxd` files, applies the
   binary masks, detects evoked events, refines nearby clusters, calculates
   event statistics, creates raster plots, and saves the results.

If `processed_data.mat` already exists in an output directory, the script skips
that recording because `ops.redo_detection` defaults to `false`. Set
`ops.redo_detection = true` to regenerate existing results.

### Expected output

Each recording is saved in its own directory under `Test_data/outputs/`, with
files including:

- `processed_data.mat`: processed pixel, ROI, event-cluster, parameter, and
  statistics data.
- `denoised_im.tif`, `ImJFig1_Max_dFoF_TopHatFiltered.tif`, and
  `ImJFig2_PxMask.tif`: intermediate image products.
- `Fig1_FirstFrame.png` through `Fig7_LabelMask.png`: preprocessing,
  screening, masking, and ROI figures.
- `Fig8_SynchroniousSynapses.png`, `Fig9_StimulationResponse.png`, and
  `rasterplot.png`: event-analysis figures where generated.
- `ROI_pxMap/` and `ROI_px/`: ROI maps and pixel-level plots where generated.

Results of matching synapses across recordings includes:

- `ClusterMatchingFig1_Distribution`, `ClusterMatchingFig2_Histogram`, and
  `ClusterMatchingFig3_Overall`: matching figures generated by
  `matching_clusters` after all recordings in the output tree have been
  processed.
- `results.mat`: containing the matched cluster summaries and combined event-cluster data. Figures are saved in
  both MATLAB `.fig` format and the configured image format.

Representative outputs already present in the repository:

![First frame](./Test_data/outputs/Cell1_1/Fig1_FirstFrame.png)

![Maximum delta F over F](./Test_data/outputs/Cell1_1/Fig2_MaxDFoF.png)

![Detected pixel mask](./Test_data/outputs/Cell1_1/Fig6_PxMask.png)

![Raster plot](./Test_data/outputs/Cell1_1/rasterplot.png)

The exact number of detected pixels, ROIs, events, and clusters depends on the
 parameter settings. Expected demo counts and
summary values: **[FILL IN after a verified clean run]**.

Expected demo run time on a normal desktop computer: **[FILL IN; record CPU,
RAM, MATLAB release, and whether a parallel pool was used]**.

## 4. Instructions for use

### Run the pipeline on new data

1. Prepare one or more single-plane time-series recordings in a
   Bio-Formats-compatible format, such as `.cxd` or `.tif`.
2. For each recording processed with `ops.use_binary_mask = true`, place a
   binary mask in the same directory as the recording. The current script
   searches for a file matching:

   ```text
   MAX_Cell*_binary_*.tif
   ```

   The mask must have the same image dimensions as the recording. Set
   `ops.use_binary_mask = false` only when automatic thresholding is intended.
3. Copy the pipeline script or create a working copy, then edit the
   **Defining parameters and file paths** section:

   ```matlab
   ops.filedir = 'PATH_TO_INPUT_FOLDER';
   ops.fileformat = '.cxd';       % Change to '.tif' or another supported format
   ops.savedir = 'PATH_TO_OUTPUT_FOLDER';
   ```
4. Set the experiment and acquisition parameters to match the recording:

   ```matlab
   ops.experiment_type = "evoked"; % or "spontaneous"
   ops.first_stim = 3;              % seconds
   ops.n_stim = 15;
   ops.stim_freq = 1;               % Hz
   ops.fs = 100;                    % frames/second
   ```

   For spontaneous recordings, set `ops.experiment_type = "spontaneous"`.
   The evoked-only stimulus parameters are then unused.
5. Review preprocessing, mask, slope/SNR, baseline-drift, ROI-area, event,
   clustering, and spike-train settings in the same section. In particular,
   set `ops.redo_detection = true` when rerunning an existing output folder.
6. Run the script after adding `Scripts/` to the MATLAB path. The script
   creates output directories automatically and writes `processed_data.mat`
   plus diagnostic figures.

### Output data for downstream analysis

Load the main result with:

```matlab
result = load(fullfile('PATH_TO_OUTPUT_FOLDER', 'RECORDING_NAME', ...
		'processed_data.mat'));
```

Depending on the run, the MAT-file contains `px`, `mask`, `ind`, `ops`, raw and
processed signals, `event_cluster`, `ROI`, and `stats`. The helper scripts
`event_stats.m`, `ROI_pxMap_allTime.m`, `show_label_mask.m`,
`show_label_mask_with_text.m`, `add_event.m`, `remove_event.m`, and
`remove_ROI.m` support analysis and manual curation of saved results.

For matching clusters across recordings, the pipeline calls
`matching_clusters(ops.savepath)` after processing a directory tree. This
expects multiple compatible `processed_data.mat` files and writes matching
figures and `results.mat` at the specified output root.

### Troubleshooting

- **`bfopen` or Java errors:** verify the `Scripts/bfmatlab/` path and the
  Bio-Formats Java configuration.
- **`tps_mlspikes` not found:** verify that `Scripts/MLspike/` is on the MATLAB
  path, including its subdirectories.
- **No mask found:** check that the mask is beside the input recording and its
  filename matches `MAX_Cell*_binary_*.tif`, or disable binary-mask use.
- **Existing data is skipped:** set `ops.redo_detection = true`.
- **Out-of-memory or slow processing:** process fewer recordings at once,
  reduce the image/time-series size, or record the available-memory guidance
  here: **[FILL IN]**.
