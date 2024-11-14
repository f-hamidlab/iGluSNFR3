%% This script detects synaptic activities from a .cxd or .tif file
% Works for:
% - single plane image 
% - iGluSNFR3 probe
%
% Dependencies:
% 1) MLspike toolbox (https://github.com/MLspike)
% 2) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html, https://bio-formats.readthedocs.io/en/v7.0.0/developers/matlab-dev.html)
% 3) Additional scripts:
%       - cluster_events.m
%       - event_stats.m
%       - findpeaks_2d.m
%       - loadBioformats.m
%       - remove_ROI.m
%       - remv_zero_padding.m
%       - save_data.m
%       - show_label_mask.m
%       - sliding_window_filter.m
%       - zero_padding.m
%
% Additional scripts/functions useful for manual curation of data:
%       - add_event.m
%       - event_stats.m
%       - plot_signal_pixel.m
%       - remove_event.m
%       - remove_ROI.m
%       - ROI_pxMap_allTime.m
%       - save_data.m
%       - show_label_mask.m
%       - spike_train.m
%       - spike_train_par.m
%
% Last modified: 2024-10-18 17:26 pm
%                (added tophat filter for stack)

close all
clear
clc

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

% path to data folder
foldername = '/PATH/TO/PROJECT/originals/';

% path to other required functions
addpath('./')
addpath('../bfmatlab/')
addpath('../MLspike/brick/')
addpath('../MLspike/spikes/')

loop_through_folder(foldername)


%%
function loop_through_folder(foldername)
    %% loops through folder and subfolders
    filelist = dir(foldername);
    
    for i = 3:length(filelist) % the first two are '.' and '..', skip
        if filelist(i).isdir % is a folder
            loop_through_folder(fullfile(filelist(i).folder,filelist(i).name));
        
        elseif contains(filelist(i).name, '.cxd') % is a .cxd
            multi_cell_activity_detection(filelist(i).folder,filelist(i).name); % processing for each file

        else % other file type
            continue
        end
    end
end

function multi_cell_activity_detection(foldername, filename)
%% Defining parameters (check every time)
% path to data
% can take .cxd or .tif as input
ops.filename = fullfile(foldername, filename);

% path to saving directory
filename = strsplit(filename,'.');
ops.savedir = fullfile(foldername, filename{1}); % folder

% additional output format for figures
ops.fig_format = '.png'; %'.tif'

% sliding window filter
ops.baseline_percentage = 0.2; % [%]
ops.sl_window = 200; % [frame]
ops.sl_window_ST = 500; % [frame] for spike train

% threshold for slope, below which the variation in intensity is ignored
ops.m_mode = 'SD'; % 'SD' or 'absolute'
ops.m_thres = 3; % or 20 for 'absolute'
ops.SNR_thres = 5;

% threshold for removing signals with lot of curvatures in baseline
ops.BL_drift_thres = 70;

% thresholds for peak detection
ops.rising_time_thres = 5; % [frames]
ops.maxISI = 8; % [frames]
ops.findpeak_window = 7; % [frames]
ops.dfof_MinPeakHeight = 4; % for each pixel, multiple of SD
ops.MinPeakWidth = 2; % for 1d findpeaks, [sampling pts]
ops.MinPeakHeight = 150; % for 2d findpeaks

% frame rate
ops.fs = 100; % [Hz]

% filter for ROI area
ops.Area_thres = 1; % [px]

% cutoff for clustering
ops.cutoff = 5; % [pixels]

% parameters for spike train 
ops = spike_train_par(ops);

%% Turn off unnecessary warnings
warning('off','signal:findpeaks:largeMinPeakHeight')

%% Load data and segmentation
tic;
disp('Loading data ...')

% create folders for saving data/figures
if ~exist(ops.savedir, 'dir')
    mkdir(ops.savedir)
end
ops.savedir_ROIpxMap = fullfile(ops.savedir,'ROI_pxMap');
if ~exist(ops.savedir_ROIpxMap, 'dir')
    mkdir(ops.savedir_ROIpxMap)
end
ops.savedir_ROIpx = fullfile(ops.savedir,'ROI_px');
if ~exist(ops.savedir_ROIpx, 'dir')
    mkdir(ops.savedir_ROIpx)
end
ops.savedir_ROIsignal = fullfile(ops.savedir,'ROI_signal');
if ~exist(ops.savedir_ROIsignal, 'dir')
    mkdir(ops.savedir_ROIsignal)
end

[im_data, ops] = loadBioFormats(ops);

im_data = double(im_data);
ops.t = ((0: 1:ops.Nt-1)/ops.fs)';

% replicate second frame to first frame to remove artifact
im_data(:,:,1) = im_data(:,:,2);
data1 = reshape(im_data, [], ops.Nt);

toc

% plot first frame
figure()
imagesc(squeeze(im_data(:,:,1)))
title(sprintf('Frame %d',1))
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
colorbar;
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig1_FirstFrame.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig1_FirstFrame',ops.fig_format))
close(gcf)

% tophat filter (background subtraction)
tic;
disp('Applying Tophat filter...')
se = strel("cuboid",[24 24 20]);
im_data = imtophat(im_data,se);
toc;

% start timer
tic;
disp('Screening pixels...')

% assuming prominent signal with little drifting
max_data = max(im_data,[],3);
median_data = median(im_data,3);
max_df = max_data-median_data;
max_dff = max_df./median_data;

figure;
imagesc(max_dff)
title('Max. delta F over F')
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig2_MaxDFoF.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig2_MaxDFoF',ops.fig_format))
close(gcf)

% tophat filtering to correct for uneven illumination
I = max_dff;
se = strel('disk',12);
tophatFiltered = imtophat(I,se);
figure
imagesc(tophatFiltered)
title('max dff tophat filtered')
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig3_TopHatFilter.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig3_TopHatFilter',ops.fig_format))
close(gcf)

% thresholding
I = tophatFiltered;
T = graythresh(I);
max_dff_BW = imbinarize(I,T);
figure;
imagesc(max_dff_BW)
title('max dff segmented')
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig4_Thresholding.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig4_Thresholding',ops.fig_format))
close(gcf)

toc

%%
tic;
disp('Shortlisting signals for analysis...')

ind = find(max_dff_BW);
idx = num2cell(ind);
px = repmat(struct('idx',1, 'df', [], 'dff_t', [], 'dff', [], 'ipt', []), length(idx), 1 );
[px.idx] =idx{:};
signal_raw = data1(ind,:)';

toc
disp('Filtering pixels with small SNR...')

% get rising slopes
[~,signal_edge] = gradient(signal_raw);

% calculate SNR
signal_edge_signal = max(signal_edge, [], 1);
signal_edge_noise = std(signal_edge);
signal_edge_SNR = signal_edge_signal./signal_edge_noise;
y = num2cell(double(signal_edge_SNR));
[px.SNR] =y{:};
clear signal_edge_noise signal_edge_signal signal_edge_SNR signal_edge

% remove ROI with SNR > threshold
k = (vertcat(px.SNR) > ops.SNR_thres);
px = px(k,:);
signal_raw = signal_raw(:,k);
ind = ind(k);

toc

% Getting baseline 
disp('Calculating baselines...')
signal_baseline = sliding_window_filter(signal_raw, ops.baseline_percentage, ops.sl_window);
% signal_baseline_movemean = movmean(signal_raw, 10, 1);% window of 10
signal_df = signal_raw - signal_baseline;
% signal_df_movemean = signal_raw - signal_baseline_movemean;
signal_baseline = signal_baseline + median(signal_df,1);
signal_df = signal_df - median(signal_df,1);
signal_dfof = signal_df./signal_baseline;
signal_dfof_movemean = movmean(signal_dfof, 5, 1); % window of 5

% Remove px with large change in baseline*******************************%
disp('Removing pixels with large baseline drift...')
polyfit_error = zeros(1,size(signal_baseline,2));
for i = 1:size(signal_baseline,2)
    y = signal_baseline(:,i);
    p = polyfit(ops.t, y, 1);
    y1 = polyval(p,ops.t);
    polyfit_error(i) = rms(y-y1);
end

figure;
plot(polyfit_error)
hold on
yline(ops.BL_drift_thres,'k','threshold')
title('Baseline Drift Threshold')
xlabel('pixel')
ylabel('Error of Linear Fit')
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig5_BaselineDriftThreshold.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig5_BaselineDriftThreshold',ops.fig_format))
% close(gcf)

% remove ROI with Error of Linear Fit > threshold
k = find(polyfit_error < ops.BL_drift_thres);
remove_px(k)

toc

% find timing of events
disp('Screening for events...')
signal_edge = Gaussian_Derivative_Filter_padded(signal_dfof, ops.m_mode, ops.m_thres);
signal_edge_nve = Gaussian_Derivative_Filter_padded(-signal_dfof, ops.m_mode, ops.m_thres);

% use parallel computing
parfor i = 1:size(signal_raw,2) % for each pixel
    disp(i)
    [pks,locs] = findpeaks(signal_edge(:,i));
    [pks_nve,locs_nve] = findpeaks(signal_edge_nve(:,i));
    % skip the first and last 5 frames
    pks = pks(locs>5);
    locs = locs(locs>5);
    pks = pks(locs<ops.Nt-5);
    locs = locs(locs<ops.Nt-5);

    px(i).slope_pk = pks;
    px(i).slope_t = locs;

    signal_dfof_std = signal_dfof(:,i);

    for j=1:length(locs)
        y=locs(j);
        signal_dfof_std(y:y+40) = NaN;

    end
    px(i).STD = std(signal_dfof_std(:),'omitnan');

    for j = 1:length(locs) % for each frame where peak slope is located
        t_start = locs(j);
        t_end = min([ops.Nt, t_start+ops.findpeak_window]);

        [~, dff_t] = findpeaks(signal_dfof(t_start:t_end,i),'NPeaks',1, ...
            'MinPeakHeight', ops.dfof_MinPeakHeight*px(i).STD, ...
            'WidthReference','halfheight',...
            'MinPeakWidth',ops.MinPeakWidth);
        
        % find the first zero-crossing
        % if the first zero-crossing is over threshold, remove pt
        if ~isempty(dff_t)
            dff_t = dff_t+t_start-1;
            
            if ~isempty(px(i).dff_t)
                last_pk = px(i).dff_t(end);
                dis = dff_t-last_pk;
                if dis<ops.maxISI
                    px(i).df(end+1) = signal_df(dff_t,i);
                    px(i).dff_t(end+1) = dff_t;
                    px(i).dff(end+1) = signal_dfof(dff_t,i);
                    continue
                else 
                    % check time of rising slope
                    STD_crossing = find(signal_dfof_movemean(1:dff_t,i)<3*px(i).STD,1,"last");
                    t_find_ch_pt = max([STD_crossing-50, round((last_pk+dff_t)/2), 1]);
                    ipt = findchangepts(signal_dfof_movemean(t_find_ch_pt:dff_t,i),Statistic="linear", MaxNumChanges=1);
                    ipt = ipt + t_find_ch_pt;
                    if ~isempty(ipt)
                        px(i).ipt(end+1) = ipt;
                    end
                    % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                    rising_time = dff_t-ipt;
                end
            else 
                % check time of rising slope
                STD_crossing = find(signal_dfof_movemean(1:dff_t,i)<3*px(i).STD,1,"last");
                ipt = findchangepts(signal_dfof_movemean(max([STD_crossing-50, 1]):dff_t,i),Statistic="linear", MaxNumChanges=1);
                ipt = ipt + max([STD_crossing-50, 1])-1;
                if ~isempty(ipt)
                    px(i).ipt(end+1) = ipt;
                end
                % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                rising_time = dff_t-ipt;
            end
            
            if rising_time>ops.rising_time_thres
                continue
            end

            px(i).df(end+1) = signal_df(dff_t,i);
            px(i).dff_t(end+1) = dff_t;
            px(i).dff(end+1) = signal_dfof(dff_t,i);

        end
    end
    px(i).n_spikes = length(px(i).dff_t);
end
toc

%% remove px with no spikes
tic;
disp('Removing pixels with no event...')
n_spikes = vertcat(px.n_spikes);
ind = vertcat(px.idx);
k = find(n_spikes>0);

pxmask = zeros(ops.Ny, ops.Nx);
pxmask(ind(k)) = 1;

figure
imagesc(pxmask)
title('max dff segmented')
xlabel('X [px]')
ylabel('Y [px]')
colormap gray
axis image
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig6_PxMask.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig6_PxMask',ops.fig_format))
close(gcf)

remove_px(k)
toc

%% Defining ROI

pxmask = zeros(ops.Ny, ops.Nx);
pxmask(ind) = 1;

tic;
disp('Defining ROI...')
L = bwlabel(pxmask,8);
ROI = regionprops(L,"PixelIdxList", "PixelList", "Area", "Centroid");

% remove ROI with Area> threshold
ROI = ROI([ROI.Area] > ops.Area_thres);

figure;
imagesc(L)
title('label mask')
xlabel('X [px]')
ylabel('Y [px]')
axis image
cmap = parula;
cmap(1,:)=[0,0,0];
colormap(cmap)
% save figure
saveas(gcf, strcat(ops.savedir,filesep, 'Fig7_LabelMask.fig'))
saveas(gcf, strcat(ops.savedir,filesep, 'Fig7_LabelMask',ops.fig_format))
close(gcf)

% initialize empty array
empty_cell = cell(length(ROI),1);
[ROI.df] = empty_cell{:};
[ROI.dfof] = empty_cell{:};
[ROI.dfft] = empty_cell{:};
[ROI.t] = empty_cell{:}; % non-repeating timing of events

empty_cell = cell(length(px),1);
[px.label] = empty_cell{:};

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
figure
% for each ROI
for n = 1:size(ROI,1)
    event_cluster = cluster_events(ROI,event_cluster,n);
end
event_cluster = event_cluster'; % correct dimension
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
                    if all(abs(t1-t2)<=1) % falls within one frame before and after
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
        if length(ROI_0) > 1 % grouping different ROIs
            ROI(end+1).Area = sum(vertcat(ROI(ROI_0).Area));
            ROI(end).PixelIdxList = vertcat(ROI(ROI_0).PixelIdxList);
            ROI(end).PixelList = vertcat(ROI(ROI_0).PixelList);
            ROI(end).Centroid = mean(ROI(end).PixelList);
            ROI = ROI_properties(ROI, length(ROI));
            % 
            % % relabel px.label
            % ia = find(ismember(label, ROI_0));
            % for k = 1:length(ia)
            %     px(ia(k)).label = length(ROI);
            % end
    
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

%% plot signal for each ROI (overall)
tic;
disp('Plotting signals for each ROI (overall)...')
figure;
set(gcf,'Position',[300 300 1500 400])
for n = vertcat(event_cluster.ROI)'
    title(sprintf('ROI: %03d, Area: %d, n spikes: %d',n, ROI(n).Area, length(ROI(n).t)))
    hold on
    plot(ops.t, ROI(n).df)
    scatter(ops.t(ROI(n).t), ROI(n).df(ROI(n).t))
    ylabel('df')
    xlabel('Time [s]')
    % save figure
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),'.fig'))
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),ops.fig_format))
    clf(gcf)

end
close(gcf)
toc

%% plot signal for each ROI (all pixels)
tic;
disp('Plotting signal for each ROI (all pixels)...')
label = vertcat(px.label);
figure;

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
            saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_%d',n,Q),'.fig'))
            saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_%d',n,Q),ops.fig_format))
            clf(gcf)
        end
    end
    
end

close(gcf)
toc

%% event statistics
stats = event_stats(event_cluster, ops);

%% save data
save_data()

%%
show_label_mask_with_text(event_cluster, ROI, ops)

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
        [~, t] = findpeaks(N_w,'MinPeakHeight',1.5);
        % for each time point, find the temporal peak in ROI(n).df
        % it could be w points before/after t
        w = 2;
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
        evalin('caller','signal_edge_nve = signal_edge_nve(:,k);');
    end
end



%% Gaussian filter
function edge = Gaussian_Derivative_Filter_padded(data, mode, thres)
    tic;
    disp('Applying 1st order Gaussian filter...')

    % sigma for guassian filter
    sigma_y = 2;
    
    % length of zero-padding
    zeropad_len = 50;
    
    % zero-padding
    edge = zero_padding(data, zeropad_len);
    
    % 1st Order Derivative 1D Gaussian filter, detects slopes in time dimension
    edge = Gaussian_Derivative_Filter(edge, sigma_y);
    
    % remove zero-padding
    edge = remv_zero_padding(edge, zeropad_len);
    
    % correct for direction of the slope
    edge = -edge;
    
    % thresholding to remove slow flucturations in signal
    switch mode
        case 'absolute'
            edge(edge<thres)=0;
        case 'SD'
            STD = std(edge);
            edge(edge<thres*STD)=0;
    end

    toc
end

function edge = Gaussian_Derivative_Filter_padded_2nd(data, mode, thres)
    tic;
    disp('Applying 1st order Gaussian filter...')

    % sigma for guassian filter
    sigma_y = 2;
    
    % length of zero-padding
    zeropad_len = 50;
    
    % zero-padding
    edge = zero_padding(data, zeropad_len);
    
    % 1st Order Derivative 1D Gaussian filter, detects slopes in time dimension
    edge = Gaussian_Derivative_Filter(edge, sigma_y);
    
    % remove zero-padding
    edge = remv_zero_padding(edge, zeropad_len);
    
    % correct for direction of the slope
    edge = -edge;
    % 
    % % thresholding to remove slow flucturations in signal
    % switch mode
    %     case 'absolute'
    %         edge(edge<thres)=0;
    %     case 'SD'
    %         STD = std(edge);
    %         edge(edge<thres*STD)=0;
    % end

    toc
end

function im = Gaussian_Derivative_Filter(im, sigma_y)
    y_values = -ceil(4*sigma_y):ceil(4*sigma_y);
    filter_y = calc1stOrderDerivative1DGaussian(y_values, sigma_y)';
    im = imfilter(im, filter_y,'symmetric');
end

function [ D_values ] = calc1stOrderDerivative1DGaussian(x_values, sigma)
    % Getting the 1D Gaussian values
    G_values = calc1DGaussian(x_values, sigma);
    
    % Obtaining the 1st order derivative of the 
    D_values = -2.*x_values./2./sigma^2.*G_values;
end

function [ G_values ] = calc1DGaussian(x_values, sigma)
    G_values = 1/sigma/sqrt(2*pi).*exp(-x_values.^2./2./sigma^2);
end


