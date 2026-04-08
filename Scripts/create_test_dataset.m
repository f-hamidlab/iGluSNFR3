%% create_test_dataset.m
% Creates synthetic test dataset by manipulating original microscopy images
%
% DESCRIPTION:
%   This script generates three test datasets from a single source image:
%   1. Original image saved as multipage TIFF
%   2. Constant frame (first frame replicated throughout entire stack)
%   3. Shifted image (original shifted by N frames with first N frames filled by first frame)
%
% DEPENDENCIES:
%   1) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html)
%   2) loadBioFormats.m (in Scripts folder)
%
% USAGE:
%   1. Set FRAME_SHIFT parameter below (default: 50 frames)
%   2. Update source and output paths as needed
%   3. Run script: generates 3 TIFF files with metadata preserved
%
% OUTPUT:
%   - Cell1_1.tif: Copy of original image
%   - Cell1_2.tif: Constant image (first frame repeated)
%   - Cell1_3.tif: Shifted image (shifted by FRAME_SHIFT frames)
%
% Last modified: 2026-04-01

clear
close all
clc

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% USER CONFIGURATION ============================================================

% ========== FRAME SHIFT PARAMETER ==========
FRAME_SHIFT = 50;  % Number of frames to shift the image (user-configurable)

% ========== FILE PATHS ==========
% Source image path (CXD format)
source_image_path = '../../../originals/iGluSNFR3 evoked 250312 Halo C/Image1/Cell1_1.cxd';

% Second source image path (CXD format) - optional
source_image_path_2 = '../../../originals/iGluSNFR3 evoked 250312 Halo C/Image1/Cell1_2.cxd';

% Output directory for test images
output_directory = '../../../originals/iGluSNFR3 evoked 250312 Halo C/Image1_test_data/';

% ========== OUTPUT FILENAMES ==========
output_filename_1 = 'Cell1_1.tif';  % Original image
output_filename_2 = 'Cell1_2.tif';  % Constant frame image
output_filename_3 = 'Cell1_3.tif';  % Shifted image
output_filename_4 = 'Cell1_4.tif';  % Second source image

% ========== ADD DEPENDENCIES TO PATH ==========
addpath('./Scripts/')                       % Core analysis scripts
addpath('./Scripts/bfmatlab/')              % Bio-Format Toolbox

%% STAGE 1: VALIDATE INPUT ======================================================

disp('Validating input parameters...')
disp(['Current working directory: ' pwd])

% Get script directory for resolving relative paths
script_dir = fileparts(mfilename('fullpath'));
disp(['Script directory: ' script_dir])

% Add Bio-Formats toolbox to path using absolute paths
addpath(fullfile(script_dir, 'bfmatlab'));
disp(['Added to path: ' fullfile(script_dir, 'bfmatlab')])

% Convert relative paths to absolute if needed
if ~isfile(source_image_path)
    abs_path_1 = fullfile(script_dir, source_image_path);
    if isfile(abs_path_1)
        source_image_path = abs_path_1;
        disp(['Resolved source_image_path to: ' source_image_path])
    else
        error('Source image file not found at:\n  Relative: %s\n  Absolute: %s', source_image_path, abs_path_1)
    end
else
    disp(['Found source image at: ' source_image_path])
end

% Check if second source image exists
if ~isfile(source_image_path_2)
    abs_path_2 = fullfile(script_dir, source_image_path_2);
    if isfile(abs_path_2)
        source_image_path_2 = abs_path_2;
        disp(['Resolved source_image_path_2 to: ' source_image_path_2])
    else
        error('Second source image file not found at:\n  Relative: %s\n  Absolute: %s', source_image_path_2, abs_path_2)
    end
else
    disp(['Found second source image at: ' source_image_path_2])
end

% Create output directory if it doesn't exist
if ~isfolder(output_directory)
    abs_output_dir = fullfile(script_dir, output_directory);
    output_directory = abs_output_dir;
    disp(['Resolved output_directory to: ' output_directory])
end

if ~isfolder(output_directory)
    mkdir(output_directory)
    disp(['Created output directory: ' output_directory])
end

% Validate FRAME_SHIFT parameter
if ~isnumeric(FRAME_SHIFT) || FRAME_SHIFT < 0 || FRAME_SHIFT ~= floor(FRAME_SHIFT)
    error('FRAME_SHIFT must be a non-negative integer')
end

disp(['Frame shift parameter: ' num2str(FRAME_SHIFT) ' frames'])
disp(' ')

%% STAGE 2: LOAD SOURCE IMAGE ===================================================

disp('Loading source image...')
tic;

% Create ops structure for loadBioFormats
ops = struct();
ops.filename = source_image_path;

% Load image using Bio-Formats (returns image_data and metadata in ops)
[im_data_original, ops] = loadBioFormats(ops, verbose=true);

% Convert to double precision for processing
im_data_original = double(im_data_original);

% Get image dimensions [height, width, frames]
[Ny, Nx, Nt] = size(im_data_original);
disp(['Image dimensions: ' num2str(Ny) ' x ' num2str(Nx) ' x ' num2str(Nt) ' (H x W x T)'])

% Get data type and bit depth of original data
original_class = class(im_data_original);
disp(['Original data class: ' original_class])

% Validate frame shift is not larger than stack
if FRAME_SHIFT >= Nt
    error('FRAME_SHIFT (%d) must be smaller than total number of frames (%d)', FRAME_SHIFT, Nt)
end

time_load = toc;
disp(['Data loading completed in ' num2str(time_load, '%.2f') ' seconds'])
disp(' ')

%% STAGE 3: CREATE IMAGE 1 - ORIGINAL (COPY) ===================================

disp('Creating Image 1: Copy of original image...')
tic;

im_data_1 = im_data_original;

% Prepare for TIFF writing (convert to appropriate data type)
% Determine output data type based on value range
min_val = min(im_data_1, [], 'all');
max_val = max(im_data_1, [], 'all');

if max_val <= 255
    im_data_1_write = uint8(im_data_1);
elseif max_val <= 65535
    im_data_1_write = uint16(im_data_1);
else
    im_data_1_write = uint32(im_data_1);
end

% Write multipage TIFF with metadata
output_path_1 = fullfile(output_directory, output_filename_1);
write_multipage_tiff(output_path_1, im_data_1_write, ops);

time_create_1 = toc;
disp(['  Image 1 saved to: ' output_path_1])
disp(['  Creation completed in ' num2str(time_create_1, '%.2f') ' seconds'])
disp(' ')

%% STAGE 4: CREATE IMAGE 2 - CONSTANT FRAME ====================================

disp('Creating Image 2: Constant frame image (first frame repeated)...')
tic;

% Extract first frame
first_frame = im_data_original(:, :, 1);

% Create image with first frame replicated
im_data_2 = repmat(first_frame, 1, 1, Nt);

% Prepare for TIFF writing
if max_val <= 255
    im_data_2_write = uint8(im_data_2);
elseif max_val <= 65535
    im_data_2_write = uint16(im_data_2);
else
    im_data_2_write = uint32(im_data_2);
end

% Write multipage TIFF with metadata
output_path_2 = fullfile(output_directory, output_filename_2);
write_multipage_tiff(output_path_2, im_data_2_write, ops);

time_create_2 = toc;
disp(['  Image 2 saved to: ' output_path_2])
disp(['  Creation completed in ' num2str(time_create_2, '%.2f') ' seconds'])
disp(' ')

%% STAGE 5: CREATE IMAGE 3 - SHIFTED IMAGE ====================================

disp(['Creating Image 3: Image shifted by ' num2str(FRAME_SHIFT) ' frames...'])
tic;

% Create output image with same dimensions as original
im_data_3 = zeros(Ny, Nx, Nt, 'double');

% Fill first FRAME_SHIFT frames with the first frame
im_data_3(:, :, 1:FRAME_SHIFT) = repmat(first_frame, 1, 1, FRAME_SHIFT);

% Shift the original image: place original frames starting at position FRAME_SHIFT+1
% This effectively creates a circular rotation where the first FRAME_SHIFT frames are empty
% and filled with the first frame, and the rest are shifted
im_data_3(:, :, FRAME_SHIFT+1:Nt) = im_data_original(:, :, 1:Nt-FRAME_SHIFT);

% Prepare for TIFF writing
if max_val <= 255
    im_data_3_write = uint8(im_data_3);
elseif max_val <= 65535
    im_data_3_write = uint16(im_data_3);
else
    im_data_3_write = uint32(im_data_3);
end

% Write multipage TIFF with metadata
output_path_3 = fullfile(output_directory, output_filename_3);
write_multipage_tiff(output_path_3, im_data_3_write, ops);

time_create_3 = toc;
disp(['  Image 3 saved to: ' output_path_3])
disp(['  Creation completed in ' num2str(time_create_3, '%.2f') ' seconds'])
disp(' ')

%% STAGE 6: LOAD AND SAVE SECOND IMAGE ==========================================

disp('Loading and saving second source image...')
tic;

% Create ops structure for loadBioFormats
ops_2 = struct();
ops_2.filename = source_image_path_2;

% Load second image using Bio-Formats
[im_data_4, ops_2] = loadBioFormats(ops_2, verbose=false);

% Convert to double precision for processing
im_data_4 = double(im_data_4);

% Get image dimensions
[Ny_2, Nx_2, Nt_2] = size(im_data_4);
disp(['  Second image dimensions: ' num2str(Ny_2) ' x ' num2str(Nx_2) ' x ' num2str(Nt_2) ' (H x W x T)'])

% Determine output data type based on value range
min_val_4 = min(im_data_4, [], 'all');
max_val_4 = max(im_data_4, [], 'all');

if max_val_4 <= 255
    im_data_4_write = uint8(im_data_4);
elseif max_val_4 <= 65535
    im_data_4_write = uint16(im_data_4);
else
    im_data_4_write = uint32(im_data_4);
end

% Write multipage TIFF with metadata
output_path_4 = fullfile(output_directory, output_filename_4);
write_multipage_tiff(output_path_4, im_data_4_write, ops_2);

time_create_4 = toc;
disp(['  Image 4 saved to: ' output_path_4])
disp(['  Creation completed in ' num2str(time_create_4, '%.2f') ' seconds'])
disp(' ')

%% SUMMARY =====================================================================

disp('================================================================================')
disp('Test dataset creation completed successfully!')
disp('================================================================================')
disp(' ')
disp('Output files:')
disp(['  1. ' output_filename_1 ' - Original image (copy)'])
disp(['  2. ' output_filename_2 ' - Constant frame image (first frame repeated)'])
disp(['  3. ' output_filename_3 ' - Shifted image (shifted by ' num2str(FRAME_SHIFT) ' frames)'])
disp(['  4. ' output_filename_4 ' - Second source image'])
disp(' ')
disp(['All files saved to: ' output_directory])
disp(' ')
disp(['Total processing time: ' num2str(time_load + time_create_1 + time_create_2 + time_create_3 + time_create_4, '%.2f') ' seconds'])

%% ============================================================================
%  HELPER FUNCTION: write_multipage_tiff
%  ============================================================================

function write_multipage_tiff(filepath, im_stack, ops)
% Write multipage TIFF file with metadata preserved
%
% INPUTS:
%   - filepath: full path to output TIFF file
%   - im_stack: 3D image array (Ny x Nx x Nt)
%   - ops: structure containing metadata from source image
%
% USES: saveastiff.m for robust TIFF writing

    % Validate input
    if ~isnumeric(im_stack) || ndims(im_stack) ~= 3
        error('Image stack must be a 3D numeric array')
    end
    
    % Save using saveastiff (supports multipage TIFF)
    options.overwrite = true;
    options.compress = 'no';
    options.message = false;
    
    res = saveastiff(im_stack, filepath, options);
    
    if res ~= 0
        error('Failed to save TIFF file: %s (error code: %d)', filepath, res)
    end
    
end

