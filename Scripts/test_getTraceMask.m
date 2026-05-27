close all
clear
clc

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

ops.filename = '../../originals/iGluSNFR3 evoked 250312 Halo C/Image1/Cell1_1.cxd/';  % Path to the data file
ops.filename_mask = '../../originals/iGluSNFR3 evoked 250312 Halo C/Image1/MAX_Cell1_1_binary_dil3.tif';  % Path to the mask file (can be same as data file if mask is stored in the same location)

% path to saving directory
ops.savedir = '../../outputs/test_data/'; % Output folder for results and figures

addpath('./Scripts/')                       % Core analysis scripts
addpath('./Scripts/bfmatlab/')              % Bio-Format Toolbox
addpath('./Scripts/frangi_filter_version2a/')    % Frangi filter

% ========== IMAGE PREPROCESSING OPTIONS ==========
ops.tophat_r = 2;                 % Tophat filter radius [pixels]
ops.min_area = 200;                % Minimum area threshold (adjust based on data)
ops.min_length = 10;              % Minimum lengthbox threshold (MaxFeretDiameter) for elongated structures (e.g., neurites)
% ops.min_aspect_ratio = 1.5;       % Minimum aspect ratio threshold (MaxFeretDiameter/MinFeretDiameter) for elongated structures (e.g., neurites)

%% Load and preprocess data
fig_handle = figure;
tiledlayout(fig_handle,3,4)


tic;
disp('Loading data ...')

% Load microscopy data using Bio-Format reader
% Returns: im_data (Ny × Nx × Nt), ops (with Ny, Nx, Nt populated)
[im_data, ops] = loadBioFormats(ops);

% Convert to double precision for numerical operations
im_data = double(im_data);

% plotting the mean projection of the original image for visualization
nexttile(1);
imagesc(mean(im_data, 3)); % Time projection for visualization
title('Original');
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');


tic;
disp('Pre_processing: de-noise ...')    

% Remove row-wise systematic noise (e.g., scanner oscillations)
% Divide each row by its low percentile value across time
p = 5;  % 5th percentile - robust to transient signals
rowMeans = prctile(im_data, p, 2);  % Ny × Nx × 1
im_data = im_data ./ rowMeans;  % Broadcasting: implicit expansion

% Remove column-wise systematic noise (often different top vs bottom)
% Many microscopes have gradient noise in Y direction

% Process top half
colMeans1 = prctile(im_data(1:ops.Ny/2, :, :), p, 1);
im_data(1:ops.Ny/2, :, :) = im_data(1:ops.Ny/2, :, :) ./ colMeans1;  % Broadcasting

% Process bottom half (separate normalization for different noise profile)
colMeans2 = prctile(im_data(ops.Ny/2+1:end,:,:), p, 1);
im_data(ops.Ny/2+1:end,:,:) = im_data(ops.Ny/2+1:end,:,:) ./ colMeans2;  % Broadcasting

% time projection to create a 2D image for visualization
im_data = mean(im_data, 3); % Mean projection along the time dimension

% plotting the mean projection of the de-noised image for visualization
nexttile(2);
imagesc(im_data); % Time projection for visualization
title('De-noised');
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

im_denoised = rescale(im_data); % Store the de-noised image for later use in mask generation
toc;


% tophat filter to extract thin structures in the image (e.g., neurites) and correct for uneven illumination
tic;
disp('Pre_processing: tophat filter ...')
se = strel('disk', ops.tophat_r);  % Structuring element for morphological operations
im_data = imtophat(im_data,se);

% plotting the tophat-filtered image for visualization
nexttile(3);
imagesc(im_data); % Time projection for visualization
title('Tophat-filtered');
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

toc;

I = im_data; % Use the tophat-filtered image for mask generation

%% Apply Frangi filter to enhance tubular structures (e.g., neurites) in the image
% This filter uses eigenvalues of the Hessian matrix to enhance line-like structures
tic;


disp('Applying Frangi filter to enhance tubular structures ...')
options.FrangiScaleRange = [1 6]; % Range of scales to analyze (adjust based on expected neurite thickness)
options.FrangiScaleRatio = 1;     % Step size between scales
options.BlackWhite = true;        % Set to true if neurites are brighter than background
I = FrangiFilter2D(I, options);
im_data = im_data - I; % Subtracting the Frangi filter output from the original to enhance the structures of interest (e.g., neurites) and suppress background noise

% plotting the Frangi-filtered image for visualization
nexttile(4);
imagesc(im_data); % Frangi filter output visualization
title('Frangi-filtered');
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

toc;

%% Otsu's thresholding to create a binary mask of the image
% This is a simple method to segment the image into foreground (e.g., neurites) and background based on intensity
tic;
disp('Creating binary mask using Otsu''s method ...')
im_data = rescale(im_data); % Rescale image to [0, 1] for Otsu's method
level = graythresh(im_data);  % Otsu's method to find optimal threshold
binary_mask = imbinarize(im_data, level);  % Create binary mask based on the threshold

% % plotting the binary mask for visualization
% nexttile(5);
% imagesc(binary_mask); % Binary mask visualization
% title('Binary Mask (Otsu''s Thresholding)');
% colormap(gray);
% axis image;
% axis off;
% xlabel('X (pixels)');
% ylabel('Y (pixels)');

toc;

%% Get mask of cell body and subtract from binary mask to get mask of neurites only
tic;
disp('Extracting cell body mask and subtracting from binary mask to isolate neurites ...')
% create binary mask from the denoised image using Otsu's thresholding to identify the cell body
level_cell_body = graythresh(im_denoised);  % Otsu's method to find optimal threshold for cell body
binary_mask_cell_body = imbinarize(im_denoised, level_cell_body);  % Create binary mask for cell body based on the threshold
% Assuming the cell body is the largest connected component in the binary mask
cc = bwconncomp(binary_mask_cell_body);
numPixels = cellfun(@numel, cc.PixelIdxList);
[~, idx] = max(numPixels); % Find the index of the largest connected component
cell_body_mask = false(size(binary_mask));
cell_body_mask(cc.PixelIdxList{idx}) = true; % Create binary mask for the cell body

% Subtract the cell body mask from the original binary mask to get the neurite mask
neurite_mask = binary_mask & ~cell_body_mask; % Logical AND to keep only neurites
% Plot both masks on one map with different colors.
% Red: cell body, Green: neurites, Yellow: overlap.
cell_neurite_overlay = zeros(size(binary_mask, 1), size(binary_mask, 2), 3);
cell_neurite_overlay(:, :, 1) = cell_body_mask;
cell_neurite_overlay(:, :, 2) = neurite_mask;

nexttile(5);
imagesc(cell_neurite_overlay);
title('Cell Body (Red) + Neurite (Green)');
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

toc;

%% image dilate to connect nearby structures and fill small gaps in the binary mask
% This helps to create more continuous masks for structures that may have been fragmented due to noise or low contrast
tic;
disp('Processing binary mask to connect nearby structures ...')

binary_mask = neurite_mask; % Use the neurite mask for further processing

se = strel('disk', 10);  % Structuring element for dilation (adjust size as needed)
binary_mask = imclose(binary_mask, se);  % Close the binary mask to connect nearby structures


se = strel('disk', 3);  % Structuring element for dilation (adjust size as needed)
binary_mask = imdilate(binary_mask, se);  % Dilate the binary mask to connect nearby structures

% plotting the dilated binary mask for visualization
nexttile(6);
imagesc(binary_mask); % Dilated binary mask visualization
title('Processed Neurite Mask');
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');
toc;


% Regionprops to extract properties of connected components in the binary mask
% and to filter out small or irrelevant regions based on size, shape, etc.
tic;
disp('Extracting region properties from binary mask ...')
cc = bwconncomp(binary_mask);
label_img = labelmatrix(cc);
stats = regionprops(cc, 'Area', 'MaxFeretProperties', 'Solidity');

% Build color-coded maps where each connected component is filled
% with its region property value.
area_map = zeros(size(binary_mask));
max_feret_map = zeros(size(binary_mask));
solidity_map = zeros(size(binary_mask));

area_vals = [stats.Area];
max_feret_vals = [stats.MaxFeretDiameter];
solidity_vals = [stats.Solidity];

for k = 1:numel(stats)
  region_idx = (label_img == k);
  area_map(region_idx) = area_vals(k);
  max_feret_map(region_idx) = max_feret_vals(k);
  solidity_map(region_idx) = solidity_vals(k);
end

filtered_mask = ismember(label_img, find(area_vals >= ops.min_area)); % Create filtered mask based on area threshold

% plotting the filtered binary mask for visualization
nexttile(6);
imagesc(filtered_mask); % Filtered binary mask visualization
title(['Filtered Binary Mask (Area > ' num2str(ops.min_area) ' pixels)']);
colormap(gray);
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

nexttile(7);
imagesc(area_map);
title('Region Area Map');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'Area (pixels)';
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

nexttile(8);
imagesc(max_feret_map);
title('Region MaxFeretDiameter Map');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'MaxFeretDiameter (pixels)';
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

nexttile(9);
imagesc(solidity_map);
title('Region Solidity Map');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'Solidity';
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

% % keep only the regions that are long and thin (e.g., neurites) based on length and aspect ratio
% filtered_mask = ismember(labelmatrix(bwconncomp(filtered_mask)), find([stats.MaxFeretDiameter] >= ops.min_length)); % Create filtered mask based on area threshold
% % aspect_ratios = [stats.MaxFeretDiameter] ./ [stats.MinFeretDiameter]; % Calculate aspect ratios
% % filtered_mask = ismember(labelmatrix(bwconncomp(filtered_mask)), find(aspect_ratios >= ops.min_aspect_ratio)); % Create filtered mask based on aspect ratio threshold

% % plotting the final filtered binary mask for visualization
% nexttile(8);
% imagesc(filtered_mask); % Final filtered binary mask visualization
% title(['Final Filtered Binary Mask (Area > ' num2str(ops.min_area) ' pixels, Length > ' num2str(ops.min_length) ' pixels)']);
% colormap(gray);
% axis image;
% axis off;
% xlabel('X (pixels)');
% ylabel('Y (pixels)');


toc;

% load the manual mask for comparison
manual_mask = imread(ops.filename_mask); % Load the manual mask image (assumed to be a binary image where the mask is represented by white pixels)
% calculate IoU between the final filtered mask and the manual mask
intersection = sum(sum(filtered_mask & manual_mask)); % Count of pixels where both masks overlap
union = sum(sum(filtered_mask | manual_mask)); % Count of pixels where either mask has a pixel
iou = intersection / union; % Intersection over Union (IoU) metric

% calculate Precision (True Positives / (True Positives + False Positives))
tp = intersection; % True Positives
fp = sum(sum(filtered_mask)) - intersection; % False Positives
precision = tp / (tp + fp); % Precision metric

disp(['IoU between filtered mask and manual mask: ' num2str(iou)]);
disp(['Precision between filtered mask and manual mask: ' num2str(precision)]);

% Create overlay visualization
% Red channel: filtered_mask only (False Positives)
% Green channel: manual_mask only (False Negatives)
% Yellow/White: both masks (True Positives)
overlay = zeros(ops.Ny, ops.Nx, 3);
overlay(:, :, 1) = filtered_mask; % Red channel for filtered mask
overlay(:, :, 2) = manual_mask;   % Green channel for manual mask
% Yellow where both overlap, Red for filtered only, Green for manual only

% plotting the overlay of both masks
nexttile(10);
imagesc(overlay); % Overlay visualization
title(sprintf('Mask Overlay (Red=Filtered, Green=Manual, Yellow=Both) - IoU: %.3f, Precision: %.3f', iou, precision));
axis image;
axis off;
xlabel('X (pixels)');
ylabel('Y (pixels)');

toc;

set(fig_handle,'Units','normalized','Position',[0 0 1 1]); % [0 0 width height]

