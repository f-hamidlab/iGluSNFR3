%% SPONTANEOUS ACTIVITY DETECTION PIPELINE
% Comprehensive pipeline for detecting synaptic events in single-plane widefield fluorescence imaging
% 
% DESCRIPTION:
%   This script analyzes time-lapse fluorescence microscopy data to detect and characterize
%   spontaneous synaptic activity using the iGluSNFR3 glutamate-sensing fluorescent reporter.
%   The pipeline performs automated segmentation, baseline correction, event detection, 
%   spatial clustering, and statistical analysis of detected events.
%
% SUPPORTED DATA:
%   - Single-plane, single-channel widefield fluorescence microscopy
%   - iGluSNFR3 fluorescent probe (glutamate imaging)
%   - Bio-Format compatible files (.cxd, .tif, .nd2, etc.)
%   - Frame rate: User-configurable (typically 100 Hz)
%
% PIPELINE STAGES:
%   1. Data loading & preprocessing (de-noising grid artifacts)
%   2. Background subtraction & tophat filtering (uneven illumination correction)
%   3. Pixel-level thresholding & signal extraction
%   4. SNR filtering (removes low-signal pixels)
%   5. Baseline drift removal (polyfit-based filtering)
%   6. Event detection (temporal derivative + peak finding per pixel)
%   7. ROI definition (spatial clustering of active pixels)
%   8. Event clustering (merges co-active nearby ROIs)
%   9. Spike train analysis & statistics
%   10. Visualization & data export
%
% OUTPUT DATA:
%   - processed_data.mat: Detected events, spike trains, ROI statistics
%   - Figures: Preprocessing, thresholding, masks, clustering, raster plots
%   - ROI_pxMap/: Pixel-level signal maps for each ROI
%   - ImJFig*.tif: ImageJ-compatible binary masks
%
% DEPENDENCIES:
%   1) MLspike toolbox (https://github.com/MLspike)
%   2) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html)
%   3) Additional scripts in ./Scripts/ folder
%
% PERFORMANCE OPTIMIZATIONS:
%   ✓ Distance matrix caching (clustering refinement 2-3x faster)
%   ✓ Visualization toggle (disable figures for 2-3x speedup)
%   Set ops.visualize = false for maximum speed on large datasets
%
% USAGE:
%   1. Modify ops.filedir and ops.savedir paths (lines 31-32)
%   2. Adjust processing parameters as needed (lines 37-115)
%   3. Run script: it will recursively process all files in filedir
%   4. Set ops.visualize = false to disable figure generation (faster)
%
% Last modified: 2026-02-03 12:31 am
% Current branch: main | Status: Production-ready with optimizations

close all
clear
clc

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

%% Defining parameters and file paths
% USER CONFIGURATION - Modify these settings for your dataset

% ========== FILE I/O PATHS ==========
% path to data ; The script will loop through all subfolders.
ops.filedir = '../../originals/test_data/'; % Input folder containing raw microscopy files
ops.fileformat = '.cxd'; % File format to process (.cxd, .tif, .nd2, etc.)

% path to saving directory
ops.savedir = '../../outputs/test_data/'; % Output folder for results and figures

% ========== ADD DEPENDENCIES TO PATH ==========
addpath('./Scripts/')                       % Core analysis scripts
addpath('./Scripts/bfmatlab/')              % Bio-Format Toolbox
addpath('./Scripts/MLspike/brick/')         % MLspike utilities
addpath('./Scripts/MLspike/spikes/')        % Spike detection functions
addpath('./Config/')                        % Configuration files  

ops = config_spontaneous(ops);              % Load default spontaneous configuration
% ops = config_evoked(ops);                   % Load default evoked configuration

%%
if ~exist(ops.savedir, 'dir')
    mkdir(ops.savedir)
end

%% Unit conversion: Convert time-based parameters from seconds to frames
% This allows all operations to work with frame indices instead of time values
ops.sl_window          = round(ops.sl_window * ops.fs);              % Baseline window [frames]
ops.sl_window_ST       = round(ops.sl_window_ST * ops.fs);           % Spike train window [frames]

% Convert spontaneous event detection thresholds to frames
ops.spontaneous.rising_time_thres  = round(ops.spontaneous.rising_time_thres * ops.fs);
ops.spontaneous.maxISI             = round(ops.spontaneous.maxISI * ops.fs);
ops.spontaneous.findpeak_window    = round(ops.spontaneous.findpeak_window * ops.fs);
ops.spontaneous.peakWidth          = round(ops.spontaneous.peakWidth * ops.fs);
ops.spontaneous.MinPeakWidth       = round(ops.spontaneous.MinPeakWidth * ops.fs);

% Convert evoked parameters if running evoked experiment
if ops.experiment_type == "evoked"
    ops.(ops.experiment_type).rising_time_thres  = round(ops.(ops.experiment_type).rising_time_thres * ops.fs);
    ops.(ops.experiment_type).maxISI             = round(ops.(ops.experiment_type).maxISI * ops.fs);
    ops.(ops.experiment_type).findpeak_window    = round(ops.(ops.experiment_type).findpeak_window * ops.fs);
    ops.(ops.experiment_type).peakWidth          = round(ops.(ops.experiment_type).peakWidth * ops.fs);
    ops.(ops.experiment_type).MinPeakWidth       = round(ops.(ops.experiment_type).MinPeakWidth * ops.fs);

    % Calculate stimulus timing in frames
    ops.stim_time = (0:ops.n_stim-1)/ops.stim_freq + ops.first_stim;
    ops.stim_frames = ops.stim_time * ops.fs;
    ops.stim_pk_search_range = arrayfun(@(x) x:x+ops.(ops.experiment_type).findpeak_window, ops.stim_frames, 'UniformOutput', false);
    ops.stim_pk_search_range = reshape(ops.stim_pk_search_range,[],1);
    ops.stim_pk_search_range = cell2mat(ops.stim_pk_search_range);
    ops.len_spike = ops.len_spike * ops.fs;
end

% Convert spike train parameters to frames
ops.ST.sumOfPeak_window = round(ops.ST.sumOfPeak_window * ops.fs);
ops.ST.gaussian_window  = round(ops.ST.gaussian_window * ops.fs);
ops.ST.gap_thres        = round(ops.ST.gap_thres * ops.fs);

%% Loop through all files in the input directory and process
loop_through_folder(ops.filedir, ops);

disp('Done:)')

%%
function loop_through_folder(foldername, ops)
    %% loops through folder and subfolders
    filelist = dir(foldername);
    ops.savepath = ops.savedir;

    for i = 1:length(filelist) 
        
        if strcmp(filelist(i).name,'.')
            continue

        elseif strcmp(filelist(i).name,'..') 
            continue
        
        elseif filelist(i).isdir % is a folder
            disp(filelist(i).name);
            
            % create new folder for saving data for each subfolder
            ops.savedir = fullfile(ops.savepath, filelist(i).name);

            if ~exist(ops.savedir, 'dir')
                mkdir(ops.savedir)  
            end
            
            loop_through_folder(fullfile(filelist(i).folder,filelist(i).name), ops);

        elseif contains(filelist(i).name, ops.fileformat) % change file format
            ops.filename = fullfile(filelist(i).folder, filelist(i).name);
            
            % create new folder for saving data for each image file
            [~,filename,~] = fileparts(ops.filename);
            ops.savedir = fullfile(ops.savepath, filename);
            if ~exist(ops.savedir, 'dir')
                mkdir(ops.savedir)  
            end

            % try
                tStart = tic;
                multi_cell_activity_detection(ops);
                tEnd = toc(tStart);
                fprintf('Pipeline completed for %s in %.2f seconds (%.2f minutes)\n', filename, tEnd, tEnd/60);
            % catch
            %     warning('Error. Check file.')
            % end
        end
    end
end

function multi_cell_activity_detection(ops)
%% Defining parameters (check every time)

[foldername,filename,~] = fileparts(ops.filename);
disp(filename);

%path to binary image
if ops.use_binary_mask
    ops.binary_file   = dir([foldername, filesep, 'MAX_Cell*_binary_*.tif']);
    ops.binary_file = fullfile(ops.binary_file.folder, ops.binary_file.name);
end

%% Turn off unnecessary warnings
warning('off','signal:findpeaks:largeMinPeakHeight')

%% Load data and segmentation


% create folders for saving data/figures
if ~exist(ops.savedir, 'dir')
    mkdir(ops.savedir)
end
ops.savedir_ROIpxMap = fullfile(ops.savedir,'ROI_pxMap');
if ~exist(ops.savedir_ROIpxMap, 'dir')
    mkdir(ops.savedir_ROIpxMap)
end

if ~ops.redo_detection
    % check if data is already processed
    output_file = fullfile(ops.savedir, 'processed_data.mat');
    if exist(output_file, 'file')
        disp("File already processed. Skipping...")
        return
    end
end

% ========== STAGE 1: DATA LOADING & TIME VECTOR ==========
tic;
disp('Loading data ...')

% Load microscopy data using Bio-Format reader
% Returns: im_data (Ny × Nx × Nt), ops (with Ny, Nx, Nt populated)
[im_data, ops] = loadBioFormats(ops);

% Convert to double precision for numerical operations
im_data = double(im_data);
% Create time vector in seconds: [0, dt, 2*dt, ..., (Nt-1)*dt]
ops.t = ((0: 1:ops.Nt-1)/ops.fs)';

% Remove first-frame artifact (often corrupted in some microscopy systems)
im_data(:,:,1) = im_data(:,:,2);

% ========== STAGE 2: PREPROCESSING - Remove systematic noise ==========
% This step removes row/column systematic noise common in some microscopes
if ops.pre_processing
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

    % Save denoised image for quality control
    options.overwrite = true;
    saveastiff(im_data, strcat(ops.savedir, filesep, 'denoised_im.tif'), options)

    % Visualize systematic noise patterns (optional)
    if ops.visualize
        fig_handle = figure;
        tiledlayout(2,2)

        nexttile
        imagesc(squeeze(rowMeans(:,1,:)))
        title("Systematic Noise: Row-wise Pattern")
        ylabel('Row Index')
        xlabel('Frame')

        nexttile
        imagesc(squeeze(colMeans1(1,:,:)))
        title("Systematic Noise: Column-wise Pattern (Top Half)")
        ylabel('Column Index')
        xlabel('Frame')

        nexttile
        imagesc(squeeze(colMeans2(1,:,:)))
        title("Systematic Noise: Column-wise Pattern (Bottom Half)")
        ylabel('Column Index')
        xlabel('Frame')

        fig_name = 'QCfig01_pre_processing';
        save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
    end

    toc;
end

% ========== Reshape data for pixel-level processing ==========
% Reshape: (Ny × Nx × Nt) → (Ny*Nx × Nt) = (n_pixels × n_timepoints)
data1 = reshape(im_data, [], ops.Nt);

toc

% ========== STAGE 3: INITIAL IMAGE ANALYSIS ==========
% Display first frame for quality inspection
if ops.visualize
    fig_handle = figure;
    imagesc(squeeze(im_data(:,:,1)))
    title(sprintf('Frame %d - Raw Intensity', 1))
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap gray
    axis image
    colorbar;    
    % save figure
    fig_name = 'Fig1_FirstFrame';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
end


if ops.bkg_subtraction
    % tophat filter (background subtraction)
    tic;
    disp('Getting background ...')
    se = strel("disk",100);
    im_data_1stFrame = imtophat(squeeze(im_data(:,:,1)),se);
    im_data_bkg = squeeze(im_data(:,:,1))-im_data_1stFrame;
    toc;
end

% start timer
tic;
disp('Screening pixels...')

% assuming prominent signal with little drifting
max_data = max(im_data,[],3);
median_data = median(im_data,3);
max_df = max_data-median_data;
max_dff = max_df./median_data;

if ops.visualize
    fig_handle = figure;
    imagesc(max_dff)
    title('Max. delta F over F')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap gray
    axis image
    % save figure
    fig_name = 'Fig2_MaxDFoF';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
end


if ops.tophat_max_dff
    % ========== STAGE 3: TOPHAT FILTER FOR UNEVEN ILLUMINATION ==========
    % Apply morphological tophat filter to correct for uneven illumination
    % (common in widefield microscopy due to vignetting)
    I = max_dff;
    se = strel('disk',ops.tophat_max_dff_r);
    tophatFiltered = imtophat(I,se);
    if ops.visualize
        fig_handle = figure;
        imagesc(tophatFiltered)
        title('Max ΔF/F: Tophat Filtered (uneven illumination corrected)')
        xlabel('X [px]')
        ylabel('Y [px]')
        colormap gray
        axis image
        
        % save figure
        fig_name = 'Fig3_TopHatFilter';
        save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
    end

    % Save as ImageJ-compatible TIFF (optional)
    if ops.save_filtered_tif
        options.overwrite = true;
        saveastiff(tophatFiltered, strcat(ops.savedir,filesep, 'ImJFig1_Max_dFoF_TopHatFiltered.tif'), options)
    end
end

if ops.use_binary_mask
    % load and apply binary mask
    BW = imread(ops.binary_file);
    BW = BW>0; % convert to binary

    if ops.tophat_max_dff
        I = tophatFiltered;
        I = imextendedmax(I, graythresh(I));
        I(~BW) = 0;
    else
        I = max_dff;
        I(~BW) = 0;
    end

    fig_handle = figure;
    imagesc(I)
    title('Applied binary mask')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap gray
    axis image
    % save figure
    fig_name = 'Fig4_BinaryMasked';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);


    % get signal overall for area within binary mask
    mask = struct('idx', 1, 'df', [], 'dff_t', [], 'dff', [], 'ipt', [], 'dff_from_ipt', []);
    mask.BW = BW;
    mask.signal_raw = mean(data1(BW,:),1)';

    if ops.experiment_type == "spontaneous"
        [mask.signal_baseline, mask.signal_df, mask.signal_dfof, mask.signal_dfof_movemean, mask.signal_edge] = signal_processing(mask.signal_raw, ops, ops.experiment_type);
        mask = px_properties(mask, ops, mask.signal_edge, mask.signal_df, mask.signal_dfof, mask.signal_dfof_movemean);
    elseif ops.experiment_type == "evoked"
        [mask.signal_baseline, mask.signal_df, mask.signal_dfof, mask.signal_dfof_movemean, mask.signal_edge] = signal_processing(mask.signal_raw, ops, "evoked_mask");
        mask = mask_properties_evoked_overall(mask, ops);
    end
    
    % plotting
    fig_handle = figure;
    sgtitle("Signal in binary mask")
    
    tiledlayout(2,1)
    nexttile
    hold on
    plot(ops.t, mask.signal_raw, 'b')
    plot(ops.t, mask.signal_baseline, 'k')
    ylabel('Intensity [A.U.]')
    xlabel('Time [s]')
    legend("raw","baseline")

    nexttile
    hold on
    plot(ops.t, mask.signal_dfof, 'b')
    scatter(ops.t(mask.dff_t), mask.dff)
    scatter(ops.t(mask.ipt), mask.signal_dfof(mask.ipt))
    yline(0)
    ylabel('\DeltaF/F')
    xlabel('Time [s]')

    % save figure
    fig_name = 'Fig4_1_MaskedDFoF';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
    
else
    mask = "Binary mask not given";
    % thresholding
    if ops.tophat_max_dff
        I = tophatFiltered;
        I = imextendedmax(I, graythresh(I));
    else
        I = max_dff;
        T = graythresh(I);
        I = imbinarize(I,T);
    end
        
    fig_handle = figure;
    imagesc(I)
    title('max dff segmented')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap gray
    axis image
    % save figure
    fig_name = 'Fig4_Thresholding';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
    

end

toc

%%
tic;
disp('Shortlisting signals for analysis...')

ind = find(I);
idx = num2cell(ind);
% pre-define fields of struct
px = repmat(struct('idx', 1, ...
                   'df', [], ...
                   'dff_t', [], ...
                   'dff', [], ...
                   'ipt', [], ...
                   'slope_pk', [], ...
                   'slope_t', [], ...
                   'STD', [], ...
                   'n_spikes', []), length(idx), 1 );
[px.idx] = idx{:};
signal_raw = data1(ind,:)';

if ops.bkg_subtraction
    % background subtraction (using implicit broadcasting)
    signal_raw = signal_raw - im_data_bkg(ind)';  % Broadcasting handles repetition
end

toc 

if ops.filter_by_slope
    % ========== STAGE 4: SNR-BASED PIXEL FILTERING ==========
    % Remove pixels with low signal-to-noise ratio (likely noise or inactive cells)
    % Focuses analysis on biologically active pixels only
    disp('Filtering pixels with small slope...')
    
    % Compute temporal gradient (proxy for signal changes/"slope")
    [~, signal_edge] = gradient(signal_raw);  % signal_edge: Nt × Npixels
    
    % Calculate signal-to-noise ratio per pixel (vectorized)
    signal_edge_noise = std(signal_edge);
    signal_edge_SNR = max(signal_edge, [], 1) ./ signal_edge_noise;  % SNR per pixel
    
    % Store SNR in pixel struct
    y = num2cell(double(signal_edge_SNR));
    [px.SNR] = y{:};
    clear signal_edge_noise signal_edge_SNR signal_edge
    
    % Keep only pixels with SNR > threshold
    k = (vertcat(px.SNR) > ops.SNR_thres);
    px = px(k,:);
    signal_raw = signal_raw(:,k);
    ind = ind(k);
    
    fprintf('Retained %d/%d pixels with SNR > %.1f\n', sum(k), length(k), ops.SNR_thres)
    toc
end

% ========== STAGE 5: BASELINE DRIFT REMOVAL (POLYFIT-BASED FILTERING) ==========
% Calculate baseline and remove pixels with large baseline drift
% Drift indicates motion artifacts, phototoxicity, or measurement instability
disp('Calculating baselines...')
[signal_baseline, signal_df, signal_dfof, signal_dfof_movemean, signal_edge] = signal_processing(signal_raw, ops, ops.experiment_type);

% Remove pixels with unstable baselines (indicate movement/noise artifacts)
disp('Removing pixels with large baseline drift...')

% Fit linear trend to baseline for each pixel (vectorized using matrix least squares)
% RMS error quantifies how much baseline drifts over recording
% Design matrix for linear regression: [t, 1]
t_design = [ops.t, ones(length(ops.t), 1)];

% Solve for all pixels simultaneously: coeffs = t_design \ signal_baseline
% This is much faster than looping and calling polyfit for each pixel
coeffs = t_design \ signal_baseline;  % 2 × Npixels matrix

% Compute fitted values for all pixels: baseline_fit = t_design * coeffs
baseline_fit = t_design * coeffs;  % Nt × Npixels

% Calculate RMS error per pixel
polyfit_error = rms(signal_baseline - baseline_fit, 1);  % 1 × Npixels

% Visualize baseline drift distribution (only if visualization enabled)
if ops.visualize
    fig_handle = figure;
    plot(polyfit_error)
    hold on
    yline(ops.BL_drift_thres, 'k', 'threshold')
    title('Baseline Drift Analysis')
    xlabel('Pixel Index')
    ylabel('RMS Error of Linear Fit')
    fig_name = 'Fig5_BaselineDriftThreshold';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
end

% Keep only pixels with stable baselines (low drift)
k = find(polyfit_error < ops.BL_drift_thres);
remove_px(k)

fprintf('Retained %d/%d pixels with stable baseline (drift < %.3f)', length(k), length(polyfit_error), ops.BL_drift_thres)
toc

% ========== STAGE 6: EVENT DETECTION (TEMPORAL DERIVATIVE + PEAK FINDING) ==========
% Detect spontaneous synaptic events at pixel level
% For each pixel: compute temporal derivative, threshold, and find peaks
disp('Screening for events...')

% PARALLEL PROCESSING: Each iteration analyzes one pixel
% This is the most computationally intensive stage
% Parallelization here provides significant speedup on multi-core systems
parfor i = 1:size(signal_raw, 2)  % For each pixel
    % Extract time series for this pixel
    px_edge = signal_edge(:, i);                % Temporal derivative
    px_df = signal_df(:, i);                    % ΔF (fluorescence change)
    px_dfof = signal_dfof(:, i);                % ΔF/F (normalized fluorescence)
    px_dfof_movemean = signal_dfof_movemean(:, i);  % Moving average of ΔF/F

    px(i) = px_properties(px(i), ops, px_edge, px_df, px_dfof, px_dfof_movemean);
    
end
toc

%% remove px with no spikes

n_spikes = vertcat(px.n_spikes);
ind = vertcat(px.idx);

if ops.remove_px_with_no_spikes
    tic;
    disp('Removing pixels with no event...')
    k = find(n_spikes>0);
    remove_px(k)
    toc
end

    pxmask = false(ops.Ny, ops.Nx);
    pxmask(ind) = true;
    
    if ops.visualize
        fig_handle = figure;
        imagesc(pxmask)
        title('max dff segmented')
        xlabel('X [px]')
        ylabel('Y [px]')
        colormap gray
        axis image
        % save figure
        fig_name = 'Fig6_PxMask';
        save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
    end
    
    % Save as ImageJ-compatible TIFF (optional)
    if ops.save_px_mask
        options.overwrite = true;
        saveastiff(pxmask, strcat(ops.savedir,filesep, 'ImJFig2_PxMask.tif'), options)
    end
    


%% Defining ROI
% ========== STAGE 7: ROI DEFINITION (SPATIAL CLUSTERING OF ACTIVE PIXELS) ==========
% Group active pixels into regions of interest (ROIs) using 8-connected components
% Each ROI typically represents a single synaptic bouton or small dendritic segment

pxmask = false(ops.Ny, ops.Nx);
pxmask(ind) = true;

tic;
disp('Defining ROI...')
L = bwlabel(pxmask,8); % 8-connected
ROI = regionprops(L,"PixelIdxList", "PixelList", "Area", "Centroid");

% remove ROI that are too big or too small
if isfield(ops, 'Area_thres_min')
    ROI = ROI([ROI.Area] > ops.Area_thres_min);
end
if isfield(ops, 'Area_thres_max')
    ROI = ROI([ROI.Area] < ops.Area_thres_max);
end

if ops.visualize
    fig_handle = figure;
    imagesc(L)
    title('label mask')
    xlabel('X [px]')
    ylabel('Y [px]')
    axis image
    cmap = parula;
    cmap(1,:)=[0,0,0];
    colormap(cmap)
    % save figure
    fig_name = 'Fig7_LabelMask';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);
end


% initialize empty array
empty_cell = cell(length(ROI),1);
[ROI.df] = empty_cell{:};
[ROI.dfof] = empty_cell{:};
[ROI.dfft] = empty_cell{:};
[ROI.t] = empty_cell{:}; % non-repeating timing of events

[px.label] = deal(0);

for n = 1:size(ROI,1)
    % find stats of corresponding pixels
    ROI = ROI_properties(ROI, n);
end

toc

%% Clustering
tic;
disp('Clustering events...')
event_cluster = struct([]);
n_cluster = 0;
if ops.plot_pxMap
    fig_handle = figure;
end
% for each ROI
for n = 1:size(ROI,1)
    event_cluster = cluster_events(ROI,event_cluster,n);
end
event_cluster = event_cluster'; % correct dimension

% remove cluster with no event
event_cluster = event_cluster([event_cluster.n_spikes] > 0);

if ops.plot_pxMap
    close(gcf)
end

toc

%% cluster ROIs by distance and time of spikes
% ========== STAGE 8: EVENT CLUSTERING (MERGES CO-ACTIVE NEARBY ROIs) ==========
% Merge ROIs that are spatially close and have synchronized spike timing
% Indicates separate detection of same synaptic event (reduces false fragmentation)
tic;
disp('Refining clusters...')

% Pre-compute distance matrix once (outside while loop)
D_matrix = compute_distance_matrix(event_cluster);

while true
    cluster_group = zeros(size(event_cluster));
    if size(event_cluster,1) ~= size(D_matrix,1)
        D_matrix = compute_distance_matrix(event_cluster);
    end

    for n = 1:size(event_cluster,1)
        for m = n+1:size(event_cluster,1)
            D = D_matrix(n,m);  % Use pre-computed distance (O(1) lookup)
            if D<ops.cutoff
                % check timing
                t1 = event_cluster(n).dff_t';
                t2 = event_cluster(m).dff_t';
    
                if length(t1) == length(t2) % same number of spikes
                    if all(abs(t1-t2)<=ops.timing_tolerance) % falls within n frame before and after
                        % create new group
                        if cluster_group(n) == 0 && cluster_group(m) == 0
                            grp = max(cluster_group) +1;
                            cluster_group(n) = grp;
                            cluster_group(m) = grp;
                        % add to existing group
                        elseif cluster_group(n) == 0
                            cluster_group(n) = cluster_group(m);
                        else % cluster_group(m) == 0
                            cluster_group(m) = cluster_group(n);
                        end
                    end
                    
                end
            end
        end
    end
    
    % group event clusters
    label = vertcat(px.label);
    for n = 1:max(cluster_group)
        i = find(cluster_group==n);
        
        % create new ROI
        ROI_0 = unique(vertcat(event_cluster(i).ROI)); % original ROIs
        i = find(ismember(vertcat(event_cluster(:).ROI), ROI_0));
        cluster_group(i) = n;

        if length(ROI_0) > 1 % grouping different ROIs
            ROI(end+1).Area = sum(vertcat(ROI(ROI_0).Area));
            ROI(end).PixelIdxList = vertcat(ROI(ROI_0).PixelIdxList);
            ROI(end).PixelList = vertcat(ROI(ROI_0).PixelList);
            ROI(end).Centroid = mean(ROI(end).PixelList);
            ROI = ROI_properties(ROI, length(ROI));
    
            % create new cluster
            event_cluster = cluster_events(ROI,event_cluster,length(ROI));
        end
    end

    event_cluster(cluster_group>0) = [];

    if length(event_cluster) == length(cluster_group)
        break % no more re-group
    end
end
toc;

%% Check n_cluster / n_px for each ROI
ops.ratio_thres = 0.02;
ops.area_thres = 200;

if isfield(event_cluster,'ROI')
    [n_cluster,~] = histcounts(vertcat(event_cluster.ROI),length(ROI));
else
    disp('No ROI detected. Skipping...')
    return
end

n_px = vertcat(ROI.Area);

ratio = n_cluster'./(n_px);
ratio(n_px<ops.area_thres) = 0;
ops.messy_ROI = find(ratio>ops.ratio_thres);

if ~isempty(ops.messy_ROI)
    warning(sprintf('Check ROI %03d \n',ops.messy_ROI'));
end

%% Spike train thresholding
% ========== STAGE 9: SPIKE TRAIN ANALYSIS & STATISTICS ==========
% Refine spike detection using spike train deconvolution or thresholding
% Produces binary spike trains for further analysis
ROI = ST_function(ROI, ops);
% adjust spike train detection parmaters (optional)
ops = spike_train_par(ops);

% add spike trains
for i = 1:length(ROI)
    if ROI(i).ST == 1
       spike_train(i);
    end 
end



%% plot signal for each ROI (overall)
if ops.plot_ROI_overall_signal
    tic;
    disp('Plotting signals for each ROI (overall)...')
    
    ops.savedir_ROIsignal = fullfile(ops.savedir,'ROI_signal');
    if ~exist(ops.savedir_ROIsignal, 'dir')
        mkdir(ops.savedir_ROIsignal)
    end
    
    fig_handle = figure;
    for n = vertcat(event_cluster.ROI)'
        title(sprintf('ROI: %03d, Area: %d, n spikes: %d',n, ROI(n).Area, length(ROI(n).t)))
        hold on
        plot(ops.t, ROI(n).df)
        scatter(ops.t(ROI(n).t), ROI(n).df(ROI(n).t))
        ylabel('df')
        xlabel('Time [s]')
        % save figure
        fig_name = sprintf('signal_ROI%03d',n);
        save_figure(fig_handle, fig_name, ops.savedir_ROIsignal, ops.fig_format, 0, [300 300 1500 400]);
        
    end
    close(fig_handle)
    toc
end

%% plot signal for each ROI (all pixels)
if ops.plot_ROI_px_signal
    tic;
    disp('Plotting signal for each ROI (all pixels)...')
    label = vertcat(px.label);
    
    ops.savedir_ROIpx = fullfile(ops.savedir,'ROI_px');
    if ~exist(ops.savedir_ROIpx, 'dir')
        mkdir(ops.savedir_ROIpx)
    end
    
    fig_handle = figure;
    for n = vertcat(event_cluster.ROI)'
        sgtitle(sprintf('ROI: %03d', n))
        n_px = length(find(label==n));
        n_row = 15;
        n_col = min([3, ceil(n_px/n_row)]);
        k = 1:n_row*n_col;
        k = reshape(k,n_row,n_col);
        k = k';
        idx = find(label==n);
    
        for i = 1:n_px
            Q = ceil(i/(n_row*n_col));
            R = rem(i,n_row*n_col);
            if R == 0
                R = n_row*n_col;
            end
            subplot(n_row, n_col, find(k==R))
            hold on
            j = idx(i);
            plot(ops.t, signal_dfof(:,j))
            scatter(ops.t(px(j).dff_t),px(j).dff)
            ylabel(j)
            if R == n_row*n_col || i == n_px
                set(gcf, 'Position', get(0, 'Screensize'));
                % save figure
                fig_name = sprintf('px_ROI%03d_%d',n,Q);
                save_figure(fig_handle, fig_name, ops.savedir_ROIpx, ops.fig_format, 0);
            end
        end
        
    end
    
    close(gcf)
    toc
end
%% event statistics
[stats, event_cluster] = event_stats(event_cluster, ops);

%% raster plot
rasterplot(ops,event_cluster);

%% save data
save_data()

%%
show_label_mask_with_text(event_cluster, ROI, ops)

close all

end

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%                                                                         %
%                            local functions                              %
%                                                                         %       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

function ROI = ROI_properties(ROI, n)
    
    ind = evalin('caller', 'ind');
    signal_df = evalin('caller', 'signal_df');
    signal_dfof = evalin('caller', 'signal_dfof');
    px = evalin('caller', 'px');
    ops = evalin('caller', 'ops');  

    % bin size for histcounts
    edges = (0:ops.Nt);
    
    % kernal for moving sum
    kernal = ones(5,1);

    % find stats of corresponding pixels
    [~,idx] = intersect(ind, ROI(n).PixelIdxList,'stable');
    for i = 1:size(idx,1)
        px(idx(i)).label = n;
        dff_t = px(idx(i)).dff_t;
        ROI(n).dfft(end+1:end+length(dff_t)) = dff_t;
    end
    
    ROI(n).df = mean(signal_df(:,idx), 2);
    ROI(n).dfof = mean(signal_dfof(:,idx), 2);


    if ROI(n).Area == 1 % only one pixel
        ROI(n).t = ROI(n).dfft;
    else
        % timing of event
        [N,~] = histcounts(ROI(n).dfft,edges);
        N_w = conv(N, kernal, 'same'); % moving sum
        if ROI(n).Area <= 3
            thres = 0.5;
        else
            thres = 1.5;
        end
        [~, t] = findpeaks(N_w,'MinPeakHeight',thres);
        % for each time point, find the temporal peak in ROI(n).df
        % it could be w points before/after t
        w = 2; % TODO
        for i = 1:length(t)
            [~,I] = max(ROI(n).df(t(i)-w:t(i)+w));
            t(i) = t(i)+I-w-1;
        end
        ROI(n).t = unique(t);
    end       
        
    assignin("caller","px",px);

end

function remove_px(k)
% keeps only pixels indexed in array k
    assignin('caller','k',k);

    evalin('caller','px = px(k,:);');
    evalin('caller','ind = ind(k);');

    evalin('caller','signal_raw = signal_raw(:,k);');

    evalin('caller','signal_baseline = signal_baseline(:,k);');
    % evalin('caller','signal_baseline_movemean = signal_baseline_movemean(:,k);');
    
    evalin('caller','signal_df = signal_df(:,k);');
    % evalin('caller','signal_df_movemean = signal_df_movemean(:,k);');

    evalin('caller','signal_dfof = signal_dfof(:,k);');
    evalin('caller','signal_dfof_movemean = signal_dfof_movemean(:,k);');
    
    % check if variable has been created
    W = evalin('caller','whos');
    doesExist = ismember('signal_edge',{W(:).name});
    if doesExist
        evalin('caller','signal_edge = signal_edge(:,k);');
        % evalin('caller','signal_edge_nve = signal_edge_nve(:,k);');
    end
end




%%
function save_figure(fig_handle, fig_name, savedir, fig_format, close_fig, fig_position)
    if nargin<6
        % set to full screen
        set(gcf, 'Position', get(0, 'Screensize'));
    else
        set(gcf,'Position',fig_position);
    end

    set(fig_handle,'Units','normalized','Position',[0 0 1 1]); % [0 0 width height]
    saveas(gcf, fullfile(savedir, [fig_name,'.fig']))
    saveas(gcf, fullfile(savedir, [fig_name, fig_format]))
    if close_fig
        close(gcf)
    end
end

function D_matrix = compute_distance_matrix(event_cluster)
    % Pre-compute distance matrix once (outside while loop)
    n_clusters = size(event_cluster, 1);
    D_matrix = zeros(n_clusters, n_clusters);
    for n = 1:n_clusters
        set1 = [event_cluster(n).x' event_cluster(n).y'];
        for m = n+1:n_clusters
            set2 = [event_cluster(m).x' event_cluster(m).y'];
            D_matrix(n,m) = min(pdist2(set1,set2), [], "all");
        end
    end
    D_matrix = D_matrix + D_matrix';  % Make symmetric for easier access
end