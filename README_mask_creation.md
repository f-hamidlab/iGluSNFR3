# Mask Creation and Manual Editing

This page explains the full workflow for creating masks automatically with `getTraceMask.m` and then refining them manually in the Napari mask editor.

## 1) Generate masks automatically in MATLAB

Use `getTraceMask.m` to process the raw microscopy files and create the initial mask outputs.

### What the script does

- Searches through the input folder recursively.
- Processes only files that match the configured filename pattern and numbering rule.
- Applies the preprocessing steps used in this project to generate a binary mask.
- Saves both a TIFF mask and a MATLAB data file for each processed image.

### Before you run it

Open `getTraceMask.m` and check these settings near the top of the file:

- `ops.filedir`: input folder containing the raw `.cxd` files.
- `ops.savedir`: output folder where the mask results will be written.
- `ops.fileformat`: file type to process, usually `.cxd`.
- `ops.index_multiplier`: which numbered files are processed.

Make sure the input and output paths match your computer and your dataset layout.

### How to run it

1. Open MATLAB.
2. Open `getTraceMask.m`.
3. Update the paths and options if needed.
4. Run the script.

### What you get

For each processed image, the script writes:

- `*_binary.tif`: a binary mask image.
- `*_mask_data.mat`: the MATLAB file used later for review and manual editing.

The `.mat` file includes the image data and mask variables used by the Napari editor.

## 2) Review and edit masks manually

After the automatic step, open the Napari plugin and edit any masks that need refinement.

For setup and usage instructions, see the plugin README:

- [Napari plugin README](Python/README_mask_modification.md)

That README explains how to create the environment, launch Napari, load the generated `*_mask_data.mat` files, and edit the filtered mask manually.

## Recommended workflow

1. Run `getTraceMask.m` to generate the initial masks.
2. Open the Napari plugin.
3. Load the generated `*_mask_data.mat` files.
4. Inspect the automatically generated filtered masks.
5. Edit masks manually where needed.
6. Save the updated mask files back to disk.

## Notes

- The automatic script is the starting point, not the final review step.
- The Napari editor is used for manual correction and fine tuning.
- Keep the MATLAB output folder and the Napari input folder aligned so the generated files are easy to find.
