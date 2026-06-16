close all
clear
clc

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

% ========== FILE I/O PATHS ==========
% The script loops recursively through filedir and processes files with index suffix k*n + 1.
ops.filedir = '../../originals/data20260616/';
ops.fileformat = '.cxd';
ops.base_filename_regex = ['^Cell.*_(\d+)', regexptranslate('escape', ops.fileformat), '$'];
ops.index_multiplier = 4;          % Files with index = (index_multiplier*n + 1) are processed, n starts from 0. For example, if index_multiplier=4, files with indices 1, 5, 9, ... are processed.

% Path to saving directory
ops.savedir = '../../outputs/data20260616/';

addpath('./Scripts/')                       % Core analysis scripts
addpath('./Scripts/bfmatlab/')              % Bio-Format Toolbox
addpath('./Scripts/frangi_filter_version2a/')    % Frangi filter

% ========== IMAGE PREPROCESSING OPTIONS ==========
ops.tophat_r = 2;                 % Tophat filter radius [pixels]
ops.min_area = 200;               % Minimum area threshold (adjust based on data)
ops.min_length = 10;              % Reserved for optional shape-based filtering
ops.visualize = true;            % true: show intermediate figure for each file
ops.save_figure = false;          % true: save intermediate figure in output folder

if ~exist(ops.savedir, 'dir')
  mkdir(ops.savedir)
end

loop_through_folder(ops.filedir, ops);
disp('Done.')

function loop_through_folder(foldername, ops)
    filelist = dir(foldername);
    ops.savepath = ops.savedir;

    for i = 1:length(filelist)
        if strcmp(filelist(i).name,'.')
                continue

        elseif strcmp(filelist(i).name,'..') 
            continue

        elseif filelist(i).isdir
            disp(filelist(i).name);
                
            % create new folder for saving data for each subfolder
            ops.savedir = fullfile(ops.savepath, filelist(i).name);

            if ~exist(ops.savedir, 'dir')
                mkdir(ops.savedir)  
            end        

            loop_through_folder(fullfile(filelist(i).folder, filelist(i).name), ops);
        
        elseif contains(filelist(i).name, ops.fileformat) && matches_index_pattern(filelist(i).name, ops)
            % Process file if it matches the format AND the filename pattern
            ops.filename = fullfile(filelist(i).folder, filelist(i).name);
            
            % create new folder for saving data for each image file
            [~,filename,~] = fileparts(ops.filename);
            ops.savedir = fullfile(ops.savepath, filename);
            if ~exist(ops.savedir, 'dir')
                mkdir(ops.savedir)  
            end

            process_single_file(ops);
        end
    end
end

function tf = matches_index_pattern(filename, ops)
    % Match files ending in an integer suffix and keep indices of form k*n + 1.
    tokens = regexp(filename, ops.base_filename_regex, 'tokens', 'once');

    if isempty(tokens)
        tf = false;
        return
    end

    idx = str2double(tokens{1});
    tf = ~isnan(idx) && idx >= 1 && ops.index_multiplier >= 1 && mod(idx - 1, ops.index_multiplier) == 0;
end

function process_single_file(ops)
    [~, basename, ~] = fileparts(ops.filename);
    fprintf('Processing %s\n', ops.filename)

    if ~exist(ops.savedir, 'dir')
        mkdir(ops.savedir)
    end

    fig_handle = [];
    if ops.visualize
        fig_handle = figure;
        tiledlayout(fig_handle, 3, 4)
    end

    % Load microscopy data using Bio-Format reader
    [im_data, ops] = loadBioFormats(ops);
    im_data = double(im_data);

    if ops.visualize
        nexttile(1);
        imagesc(mean(im_data, 3));
        title('Original');
        colormap(gray);
        axis image;
        axis off;
    end

    % Remove row/column systematic noise
    p = 5;
    rowMeans = prctile(im_data, p, 2);
    im_data = im_data ./ rowMeans;

    colMeans1 = prctile(im_data(1:ops.Ny/2, :, :), p, 1);
    im_data(1:ops.Ny/2, :, :) = im_data(1:ops.Ny/2, :, :) ./ colMeans1;

    colMeans2 = prctile(im_data(ops.Ny/2+1:end, :, :), p, 1);
    im_data(ops.Ny/2+1:end, :, :) = im_data(ops.Ny/2+1:end, :, :) ./ colMeans2;

    im_data = mean(im_data, 3);
    im_denoised = rescale(im_data);

    if ops.visualize
        nexttile(2);
        imagesc(im_data);
        title('De-noised');
        colormap(gray);
        axis image;
        axis off;
    end

    % Tophat filter to correct uneven illumination
    se = strel('disk', ops.tophat_r);
    im_data = imtophat(im_data, se);

    if ops.visualize
        nexttile(3);
        imagesc(im_data);
        title('Tophat-filtered');
        colormap(gray);
        axis image;
        axis off;
    end

    % Frangi enhancement for tubular structure suppression
    I = im_data;
    options.FrangiScaleRange = [1 6];
    options.FrangiScaleRatio = 1;
    options.BlackWhite = true;
    I = FrangiFilter2D(I, options);
    im_data = im_data - I;

    if ops.visualize
        nexttile(4);
        imagesc(im_data);
        title('Frangi-filtered');
        colormap(gray);
        axis image;
        axis off;
    end

    % Otsu thresholding
    im_filtered = rescale(im_data);
    level = graythresh(im_filtered);
    binary_mask = imbinarize(im_filtered, level);

    while sum(binary_mask(:)) < 1000
        level = level * 0.5;  % Adjust threshold if too few pixels are detected
        binary_mask = imbinarize(im_filtered, level);
    end

    % Estimate and subtract cell body mask
    level_cell_body = graythresh(im_denoised);
    binary_mask_cell_body = imbinarize(im_denoised, level_cell_body);

    cc_body = bwconncomp(binary_mask_cell_body);
    cell_body_mask = false(size(binary_mask));
    if cc_body.NumObjects > 0
        numPixels = cellfun(@numel, cc_body.PixelIdxList);
        [~, idx] = max(numPixels);
        cell_body_mask(cc_body.PixelIdxList{idx}) = true;
    end

    se = strel('disk', 10);
    cell_body_mask = imopen(cell_body_mask, se);
   
    neurite_mask = binary_mask & ~cell_body_mask;  
  
    if ops.visualize
        overlay = zeros(size(binary_mask, 1), size(binary_mask, 2), 3);
        overlay(:, :, 1) = cell_body_mask;
        overlay(:, :, 2) = neurite_mask;
        nexttile(5);
        imagesc(overlay);
        title('Cell Body (Red) + Neurite (Green)');
        axis image;
        axis off;
    end

    % Post-process mask morphology
    binary_mask = neurite_mask;
    se = strel('disk', 10);
    binary_mask = imclose(binary_mask, se);
    se = strel('disk', 3);
    binary_mask = imdilate(binary_mask, se);

    if ops.visualize
        nexttile(6);
        imagesc(binary_mask);
        title('Processed Neurite Mask');
        colormap(gray);
        axis image;
        axis off;
    end

    % Filter connected regions by area
    cc = bwconncomp(binary_mask);
    label_img = labelmatrix(cc);
    stats = regionprops(cc, 'Area', 'MaxFeretProperties', 'Solidity');

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

    filtered_mask = ismember(label_img, find(area_vals >= ops.min_area));

    if ops.visualize
        nexttile(7);
        imagesc(filtered_mask);
        title(['Filtered Binary Mask (Area > ' num2str(ops.min_area) ' pixels)']);
        colormap(gray);
        axis image;
        axis off;

        nexttile(8);
        imagesc(area_map);
        title('Region Area Map');
        colormap(gca, turbo);
        cb = colorbar;
        cb.Label.String = 'Area (pixels)';
        axis image;
        axis off;

        nexttile(9);
        imagesc(max_feret_map);
        title('Region MaxFeretDiameter Map');
        colormap(gca, turbo);
        cb = colorbar;
        cb.Label.String = 'MaxFeretDiameter (pixels)';
        axis image;
        axis off;

        nexttile(10);
        imagesc(solidity_map);
        title('Region Solidity Map');
        colormap(gca, turbo);
        cb = colorbar;
        cb.Label.String = 'Solidity';
        axis image;
        axis off;
    end

    % Save output as (original filename)_binary.tif
    output_mask_file = fullfile(ops.savedir, [basename '_binary.tif']);
    imwrite(uint8(filtered_mask) * 255, output_mask_file);
    fprintf('Saved: %s\n', output_mask_file)

    if ops.visualize && ops.save_figure
        set(fig_handle, 'Units', 'normalized', 'Position', [0 0 1 1]);
        saveas(fig_handle, fullfile(ops.savedir, [basename '_mask_steps.fig']))
        saveas(fig_handle, fullfile(ops.savedir, [basename '_mask_steps.png']))
    end

    if ~isempty(fig_handle) && isgraphics(fig_handle)
        close(fig_handle)
    end

    % save data in .mat file for downstream analysis
    output_mat_file = fullfile(ops.savedir, [basename '_mask_data.mat']);
    save(output_mat_file, 'binary_mask', 'filtered_mask', 'im_denoised', 'im_filtered', 'area_map', 'max_feret_map', 'solidity_map');  
end

