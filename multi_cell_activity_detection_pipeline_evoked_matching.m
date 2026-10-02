%% iGluSNFR3 ACTIVITY DETECTION PIPELINE
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
% Last modified: 2026-09-30 16:19
% Current branch: publication_branch
% Note: This branch is no longer updated. For improved functionalities, use the main branch instead.

close all
clear
clc

%% Defining parameters and file paths
% MODIFY HERE

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

% ========== IMAGE PREPROCESSING OPTIONS ==========
ops.pre_processing = true;                  % Remove systematic grid line noise (microscope artifact)
ops.bkg_subtraction = false;                % Background subtraction in time domain
ops.tophat_max_dff  = true;                 % Tophat filter for uneven illumination correction
ops.tophat_max_dff_r = 5;                   % Tophat filter radius [pixels]
ops.use_binary_mask = true;                 % Use user-provided binary mask (ImageJ)
ops.remove_px_with_no_spikes = true;        % Remove inactive pixels (set false for low SNR data)
ops.redo_detection = false;                 % Reprocess if results already exist

 % ========== VISUALIZATION & OUTPUT OPTIONS ==========
ops.plot_pxMap = false;                     % Plot pixel map for each ROI (memory intensive)
ops.plot_ROI_overall_signal = false;        % Plot averaged signal per ROI (memory intensive)
ops.plot_ROI_px_signal = false;             % Plot individual pixel signals per ROI (slow)
ops.fig_format = '.png';                    % Figure format (.png, .pdf, .fig)
ops.close_fig = true;                       % Close figures after saving (memory efficient)

% ========== TEMPORAL FILTERING PARAMETERS ==========
ops.baseline_percentage = 0.2;              % Baseline percentile for drift estimation [0-1]
ops.sl_window           = 3;                % Sliding window for baseline [s]
ops.sl_window_ST        = 5;                % Sliding window for spike train [s]

% ========== PIXEL-LEVEL FILTERING THRESHOLDS ==========
ops.filter_by_slope = true;                 % Enable SNR-based pixel filtering
ops.m_mode    = 'SD';                       % Mode: 'SD' (standard deviation) or 'absolute'
ops.m_thres   = 3;                          % Threshold multiplier (3 = 3x std above baseline)
ops.SNR_thres = 5;                          % Signal-to-noise ratio threshold for pixel inclusion

% ========== BASELINE DRIFT DETECTION ==========
ops.BL_drift_thres = 0.05;                  % Max RMS error for linear baseline fit (remove if exceeded)

% ========== SPONTANEOUS EVENT DETECTION PARAMETERS ==========
% These thresholds define event characteristics for spontaneous activity
ops.spontaneous.rising_time_thres  = 0.05;  % Minimum event rise time [s] - filters slow drifts
ops.spontaneous.maxISI             = 0.08;  % Maximum inter-spike interval [s] - fusion threshold
ops.spontaneous.findpeak_window    = 0.07;  % Time window for peak detection [s]
ops.spontaneous.dfof_MinPeakHeight = 3;     % Minimum peak height [multiples of pixel-level std]
ops.spontaneous.peakWidth          = 0.4;   % Expected event width [s] in ΔF/F signal
ops.spontaneous.MinPeakWidth       = 0.02;  % Minimum peak width for 1D detection [s]
ops.spontaneous.MinPeakHeight      = 0.01;  % Minimum peak height for 2D spatial detection

% ========== EXPERIMENT TYPE & STIMULATION PARAMETERS ==========
ops.experiment_type = "evoked";             % Type: "spontaneous" or "evoked"
% NOTE: evoked parameters below are unused for spontaneous experiments % 15AP at 1Hz starting at 3s
ops.first_stim      =  3;                   % Time of first stimulus [s]
ops.n_stim          = 15;                   % Total number of stimuli
ops.stim_freq       =  1;                   % Stimulation frequency [Hz]
ops.len_spike       =  2;                   % Expected response duration [s] from stim onset

% ========== EVOKED EVENT DETECTION PARAMETERS (if experiment_type == "evoked") ==========
ops.evoked.rising_time_thres  = 0.15;       % Minimum response rise time [s]
ops.evoked.maxISI             = 0.08;       % Maximum inter-spike interval [s]
ops.evoked.findpeak_window    = 0.15;       % Time window for peak search [s]
ops.evoked.dfof_MinPeakHeight = 2.5;        % Minimum response amplitude [multiples of std]
ops.evoked.peakWidth          = 2;          % Expected response width [s]
ops.evoked.MinPeakWidth       = 0.02;       % Minimum peak width for 1D detection [s]
ops.evoked.MinPeakHeight      = 0.1;        % Minimum peak height for 2D detection

% ========== ACQUISITION & ANALYSIS PARAMETERS ==========
ops.fs = 100;                               % Frame rate [Hz] - MUST match data acquisition rate

% ROI SIZE FILTERING (optional - comment out to disable)
ops.Area_thres_max = 50;                    % Maximum ROI area [pixels] - filters noise clusters
% ops.Area_thres_min = 1;                   % Minimum ROI area [pixels] - uncomment if needed

% ========== SPATIAL CLUSTERING PARAMETERS ==========
ops.cutoff = 5;                             % Distance cutoff for clustering [pixels]
                                            % ROIs within this distance + synchronized timing = merged

% ========== SPIKE TRAIN ANALYSIS PARAMETERS ==========
ops.ST.option           = "findpeaks";      % Method: "MLspike" (recommended) or "findpeaks"
ops.ST.MinPeakHeight    = 0.05;             % Minimum spike amplitude
ops.ST.sumOfPeak_window = 0.5;              % Window for spike summation [s]
ops.ST.gaussian_window  = 0.1;              % Gaussian smoothing window [s]
ops.ST.gap_thres        = 0.3;              % Gap threshold for spike identification [s]

%%

ops = spike_train_par(ops);

if ~exist(ops.savedir, 'dir')
    mkdir(ops.savedir)
end

%% unit conversion
ops.sl_window          = round(ops.sl_window * ops.fs); % [frame]
ops.sl_window_ST       = round(ops.sl_window_ST * ops.fs); % [frame] for spike train

ops.spontaneous.rising_time_thres  = round(ops.spontaneous.rising_time_thres * ops.fs); % [frame]
ops.spontaneous.maxISI             = round(ops.spontaneous.maxISI * ops.fs); % [frame]
ops.spontaneous.findpeak_window    = round(ops.spontaneous.findpeak_window * ops.fs); % [frame] for each pixel, findpeak in dfof
ops.spontaneous.peakWidth          = round(ops.spontaneous.peakWidth * ops.fs); % [frame] in time domain, for each pixel
ops.spontaneous.MinPeakWidth       = round(ops.spontaneous.MinPeakWidth * ops.fs); % [frame] for 1d findpeaks, width calculated as half-height of peak

if ops.experiment_type == "evoked"
    ops.(ops.experiment_type).rising_time_thres  = round(ops.(ops.experiment_type).rising_time_thres * ops.fs); % [frame]
    ops.(ops.experiment_type).maxISI             = round(ops.(ops.experiment_type).maxISI * ops.fs); % [frame]
    ops.(ops.experiment_type).findpeak_window    = round(ops.(ops.experiment_type).findpeak_window * ops.fs); % [frame] for each pixel, findpeak in dfof
    ops.(ops.experiment_type).peakWidth          = round(ops.(ops.experiment_type).peakWidth * ops.fs); % [frame] in time domain, for each pixel
    ops.(ops.experiment_type).MinPeakWidth       = round(ops.(ops.experiment_type).MinPeakWidth * ops.fs); % [frame] for 1d findpeaks, width calculated as half-height of peak


    ops.stim_time = (0:ops.n_stim-1)/ops.stim_freq + ops.first_stim;
    ops.stim_frames = ops.stim_time * ops.fs; % frame number
    ops.stim_pk_search_range = arrayfun(@(x) x:x+ops.(ops.experiment_type).findpeak_window, ops.stim_frames, 'UniformOutput', false);
    ops.stim_pk_search_range = reshape(ops.stim_pk_search_range,[],1);
    ops.stim_pk_search_range = cell2mat(ops.stim_pk_search_range);
    ops.len_spike = ops.len_spike * ops.fs; % frame number
end

ops.ST.sumOfPeak_window = round(ops.ST.sumOfPeak_window * ops.fs); % [frame]
ops.ST.gaussian_window  = round(ops.ST.gaussian_window * ops.fs); % [frame]
ops.ST.gap_thres        = round(ops.ST.gap_thres * ops.fs); % [frame]

%%
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
                multi_cell_activity_detection(ops);
            % catch
            %     warning('Error. Check file.')
            % end
        end
        
        % last file in the folder already processed &&
        % no subfolder in the dir
        if i == length(filelist)
            filelist(strcmp({filelist.name}, '.')) = [];
            filelist(strcmp({filelist.name}, '..')) = [];
            if sum([filelist.isdir]) == 0
                % check for matching event clusters
                matching_clusters(ops.savepath)
                %matching_clusters_with_image4(ops.savepath)

            end
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

tic;
disp('Loading data ...')

% load raw data
[im_data, ops] = loadBioFormats(ops);

im_data = double(im_data);
ops.t = ((0: 1:ops.Nt-1)/ops.fs)';

% replicate second frame to first frame to remove artifact
im_data(:,:,1) = im_data(:,:,2);

% Preprocessing: remove "grid line" noise from image
if ops.pre_processing
    tic;
    disp('Pre_processing: de-noise ...')    
    % for rows
    p = 5;
    rowMeans = prctile(im_data, p, 2);
    im_data = im_data ./ repmat(rowMeans, [1, ops.Nx, 1]);
    % for columns, the microscope used has diff. noise level for top and bottom
    % part of the image in each frame
    colMeans1 = prctile(im_data(1:ops.Ny/2, :, :), p, 1);
    im_data(1:ops.Ny/2, :, :) = im_data(1:ops.Ny/2, :, :)./ repmat(colMeans1, [ops.Ny/2, 1, 1]);
    colMeans2 = prctile(im_data(ops.Ny/2+1:end,:,:), p, 1);
    im_data(ops.Ny/2+1:end,:,:) = im_data(ops.Ny/2+1:end,:,:)./ repmat(colMeans2, [ops.Ny/2, 1, 1]);

    % save image as tif
    options.overwrite = true;
    saveastiff(im_data, strcat(ops.savedir, filesep, 'denoised_im.tif'), options)

    fig_handle = figure;
    tiledlayout(2,2)

    nexttile
    imagesc(squeeze(rowMeans(:,1,:)))
    title("systematic noise: row")
    ylabel('Row')
    xlabel('Frame')

    nexttile
    imagesc(squeeze(colMeans1(1,:,:)))
    title("systematic noise: column (top)")
    ylabel('Column')
    xlabel('Frame')

    nexttile
    imagesc(squeeze(colMeans2(1,:,:)))
    title("systematic noise: column (bottom)")
    ylabel('Column')
    xlabel('Frame')

    % save figure
    fig_name = 'QCfig01_pre_processing';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);

    toc;
    
end

%
data1 = reshape(im_data, [], ops.Nt);

toc

% plot first frame
fig_handle = figure;
first_frame = squeeze(im_data(:,:,1));
imagesc(first_frame)
title(sprintf('Frame %d',1))
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
colorbar;
% save figure
fig_name = 'Fig1_FirstFrame';
save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);


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


if ops.tophat_max_dff
    % tophat filtering to correct for uneven illumination
    I = max_dff;
    se = strel('disk',ops.tophat_max_dff_r);
    tophatFiltered = imtophat(I,se);
    fig_handle = figure;
    imagesc(tophatFiltered)
    title('max dff tophat filtered')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap gray
    axis image
    
    % save figure
    fig_name = 'Fig3_TopHatFilter';
    save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);

    % save as tiff
    options.overwrite = true;
    saveastiff(tophatFiltered, strcat(ops.savedir,filesep, 'ImJFig1_Max_dFoF_TopHatFiltered.tif'), options)
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
[px.idx] =idx{:};
signal_raw = data1(ind,:)';

if ops.bkg_subtraction
    % background subtraction
    signal_raw = signal_raw - repmat(im_data_bkg(ind)',ops.Nt,1);
end

toc 

if ops.filter_by_slope
    disp('Filtering pixels with small slope...')
    
    % get rising slopes
    [~,signal_edge] = gradient(signal_raw);
    
    % calculate SNR
    signal_edge_signal = max(signal_edge, [], 1);
    signal_edge_noise = std(signal_edge);
    signal_edge_SNR = signal_edge_signal./signal_edge_noise;
    y = num2cell(double(signal_edge_SNR));
    [px.SNR] =y{:};
    clear signal_edge_noise signal_edge_signal signal_edge_SNR signal_edge
    
    % remove ROI with SNR < threshold
    k = (vertcat(px.SNR) > ops.SNR_thres);
    px = px(k,:);
    signal_raw = signal_raw(:,k);
    ind = ind(k);
    
    toc
end

% Getting baseline 
disp('Calculating baselines...')
[signal_baseline, signal_df, signal_dfof, signal_dfof_movemean, signal_edge] = signal_processing(signal_raw, ops, ops.experiment_type);


% Remove px with large change in baseline*******************************%
disp('Removing pixels with large baseline drift...')
polyfit_error = zeros(1,size(signal_baseline,2));
for i = 1:size(signal_baseline,2)
    y = signal_baseline(:,i);
    p = polyfit(ops.t, y, 1);
    y1 = polyval(p,ops.t);
    polyfit_error(i) = rms(y-y1);
end

fig_handle = figure;
plot(polyfit_error)
hold on
yline(ops.BL_drift_thres,'k','threshold')
title('Baseline Drift Threshold')
xlabel('pixel')
ylabel('Error of Linear Fit')
% save figure
fig_name = 'Fig5_BaselineDriftThreshold';
save_figure(fig_handle, fig_name, ops.savedir, ops.fig_format, ops.close_fig);

% remove ROI with Error of Linear Fit > threshold
k = find(polyfit_error < ops.BL_drift_thres);
remove_px(k)

toc

% find timing of events
disp('Screening for events...')
% signal_edge = Gaussian_Derivative_Filter_padded(signal_dfof, ops.m_mode, ops.m_thres);

% use parallel computing
parfor i = 1:size(signal_raw,2) % for each pixel
    disp(i)
    
    px_edge = signal_edge(:,i);
    px_df = signal_df(:,i);
    px_dfof = signal_dfof(:,i);
    px_dfof_movemean = signal_dfof_movemean(:,i);

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

    pxmask = zeros(ops.Ny, ops.Nx);
    pxmask(ind) = 1;
    
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
    
    % save as tiff
    options.overwrite = true;
    saveastiff(pxmask, strcat(ops.savedir,filesep, 'ImJFig2_PxMask.tif'), options)
    


%% Defining ROI

pxmask = zeros(ops.Ny, ops.Nx);
pxmask(ind) = 1;

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
fig_handle = figure;
% for each ROI
for n = 1:size(ROI,1)
    event_cluster = cluster_events(ROI,event_cluster,n);
end
event_cluster = event_cluster'; % correct dimension

% remove cluster with no event %RJ commented out (FOV with no ROIs-271025)
%event_cluster = event_cluster([event_cluster.n_spikes] > 0);

close(gcf)

toc

%% cluster ROIs by distance and time of spikes
tic;
disp('Refining clusters...')
while true
    cluster_group = zeros(size(event_cluster));
    for n = 1:size(event_cluster,1)
        set1 = [event_cluster(n).x' event_cluster(n).y'];
        for m = n+1:size(event_cluster,1)
            set2 = [event_cluster(m).x' event_cluster(m).y'];
            D = min(pdist2(set1,set2),[],"all");
            if D<ops.cutoff
                % check timing
                t1 = event_cluster(n).dff_t';
                t2 = event_cluster(m).dff_t';
    
                if length(t1) == length(t2) % same number of spikes
                    if all(abs(t1-t2)<=2) % falls within two frame before and after
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

%%
labelMask = show_label_mask_with_text(event_cluster, ROI, ops);

%% save data
save_data()

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
