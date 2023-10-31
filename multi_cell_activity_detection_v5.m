%% This script detects synaptic activities from a .cxd or .tif file
% Works for:
% - single plane image 
% - iGluSNFR3 probe
%
% Dependencies:
% 1) MLspike toolbox (https://github.com/MLspike)
% 2) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html, https://bio-formats.readthedocs.io/en/v7.0.0/developers/matlab-dev.html)
% 3) Additional scripts:
%       - event_stats.m
%       - findpeaks_2d.m
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


close all
clear
clc

%% Defining parameters (check every time)
% path to data
% can take .cxd or .tif as input
ops.filepath = '/PATH/TO/PROJECT/originals/capture.cxd';

% path to saving directory
ops.savedir = '/PATH/TO/PROJECT/outputs/capture'; % folder

% path to other required functions
addpath('/PATH/TO/PROJECT/codes') % folder

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
ops.MinPeakHeight = 150; % for 2d findpeaks

% frame rate
ops.fs = 100; % [Hz]

% cutoff for clustering
ops.cutoff = 3; % [pixels]

% parameters for spike train 
ops = spike_train_par(ops);

%% Turn off unnecessary warnings
warning('off','signal:findpeaks:largeMinPeakHeight')

%% Load data and segmentation
tic;
disp('Loading data ...')

[~,v] = bfCheckJavaPath();  % added such that path to Bio-Format Toolbox is know
data = bfopen(ops.filepath);

% create folders for saving data/figures
if ~exist(ops.savedir, 'dir')
    mkdir(ops.savedir)
end
ops.savedir_ROIpxMap = strcat(ops.savedir,filesep,'ROI_pxMap');
if ~exist(ops.savedir_ROIpxMap, 'dir')
    mkdir(ops.savedir_ROIpxMap)
end
ops.savedir_ROIpx = strcat(ops.savedir,filesep,'ROI_px');
if ~exist(ops.savedir_ROIpx, 'dir')
    mkdir(ops.savedir_ROIpx)
end
ops.savedir_ROIsignal = strcat(ops.savedir,filesep,'ROI_signal');
if ~exist(ops.savedir_ROIsignal, 'dir')
    mkdir(ops.savedir_ROIsignal)
end

data = data{1,1};
data = data(:,1);

ops.Nx = size(data{1},2);
ops.Ny = size(data{1},1);
ops.Nt = size(data,1);
ops.t = ((0: 1:ops.Nt-1)/ops.fs)';

data = cell2mat(data);
data = reshape(data,ops.Ny,ops.Nt,ops.Nx);
data = permute(data, [1,3,2]);
data = double(data);

% replicate second frame to first frame to remove artifact
data(:,:,1) = data(:,:,2);
data1 = reshape(data, [], ops.Nt);

toc

% plot first frame
figure()
imagesc(squeeze(data(:,:,1)))
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

% start timer
tic;
disp('Screening pixels...')

% assuming prominent signal with little drifting
max_data = max(data,[],3);
median_data = median(data,3);
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
signal_baseline_movemean = movmean(signal_raw, 10, 1);% window of 10
signal_df = signal_raw - signal_baseline;
signal_df_movemean = signal_raw - signal_baseline_movemean;
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
            'MinPeakHeight', 3*px(i).STD);
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
tic;
disp('Defining ROI...')
L = bwlabel(pxmask,8);
ROI = regionprops(L,"PixelIdxList", "PixelList", "Area", "Centroid");

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


% bin size for histcounts
edges = (0:ops.Nt);

% kernal for moving sum
kernal = ones(5,1);

% initialize empty array
empty_cell = cell(length(ROI),1);
[ROI.df] = empty_cell{:};
[ROI.dfft] = empty_cell{:};
[ROI.t] = empty_cell{:}; % non-repeating timing of events

empty_cell = cell(length(px),1);
[px.label] = empty_cell{:};

for n = 1:size(ROI,1)
    % find stats of corresponding pixels
    [~,idx] = intersect(ind, ROI(n).PixelIdxList,'stable');
    for i = 1:size(idx,1)
        px(idx(i)).label = n;
        dff_t = px(idx(i)).dff_t;
        ROI(n).dfft(end+1:end+length(dff_t)) = dff_t;
    end
    
    ROI(n).df = mean(signal_df(:,idx), 2);
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
end

toc

%% plot signal for each ROI (overall)
tic;
disp('Plotting signals for each ROI (overall)...')
figure;
set(gcf,'Position',[300 300 1500 400])
for n = 1:size(ROI,1)
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

for n = 1:max(label)
    sgtitle(sprintf('ROI: %03d', n))
    n_px = length(find(label==n));
    n_row = 15;
    n_col = ceil(n_px/n_row);
    k = 1:n_row*n_col;
    k = reshape(k,n_row,n_col);
    k = k';
    idx = find(label==n);

    for i = 1:n_px
        subplot(n_row, n_col, find(k==i))
        hold on
        j = idx(i);
        plot(ops.t, signal_dfof(:,j))
        scatter(ops.t(px(j).dff_t),px(j).dff)
        ylabel(j)
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    % save figure
    saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d',n),'.fig'))
    saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d',n),ops.fig_format))
    clf(gcf)
end

close(gcf)
toc

%% Clustering
tic;
disp('Clustering events...')
% preallocate structure for speed
N = size(ROI,1)*3;    % assuming a maximum of 3 clusters per ROI
event_cluster = repmat(struct('ROI',[]), N, 1 );
n_cluster = 0;
figure
% for each ROI
for n = 1:size(ROI,1)
    fprintf('ROI%03d\n',n)
    
    % find stats of corresponding pixels
    [~,idx] = intersect(ind, ROI(n).PixelIdxList,'stable');

    t = ROI(n).t;
    
    % get location of each spike
    PixelList = ROI(n).PixelList;
    xLim = [min(PixelList(:,1))-5, max(PixelList(:,1))+5];
    yLim = [min(PixelList(:,2))-5, max(PixelList(:,2))+5];
    x = xLim(1):xLim(2);
    y = yLim(1):yLim(2);
    [X,Y] = meshgrid(x,y);

    df_map0 = zeros(yLim(2)-yLim(1)+1,xLim(2)-xLim(1)+1);
    T = [];

    if ROI(n).Area == 1 % only one pixel
        peak_x = PixelList(:,2);
 
    else % calculate location of event
        peak_x = [];
        peak_y = [];
        peak_amp = [];
        peak_t = [];
        for k = 1: length(t)
            df_map = df_map0;
            
            for i = 1:size(idx,1)
                x_loc = PixelList(i,1)-xLim(1)+1;
                y_loc = PixelList(i,2)-yLim(1)+1;
                df_map(y_loc, x_loc) = signal_df(t(k),idx(i));
            end
            df_map(df_map<0)=0;
            z = medfilt2(df_map);
            peak_inds = findpeaks_2d(X, z, ops.MinPeakHeight);
            
            n_peak = sum(peak_inds,'all');

            if n_peak <= 1 % if there is only one spatial peak
                signal_dfof1 = signal_dfof(t(k),idx(:));
                peak_amp(end+1) = max(signal_dfof1);
        
                signal_df1 = signal_df(t(k),idx(:));
                signal_df1(signal_df1<max(signal_df1)/2) = 0;
                
                % weighted average
                peak_x(end+1) = sum(PixelList(:,1).*signal_df1')/sum(signal_df1);
                peak_y(end+1) = sum(PixelList(:,2).*signal_df1')/sum(signal_df1);

                peak_t(end+1) = t(k);

            elseif n_peak ~= 0 % more than one peak
                peak_inds = sub2ind([ops.Ny, ops.Nx], Y(peak_inds), X(peak_inds));
                for i = 1:n_peak
                    df_map_i = df_map;
                    [row, col] = ind2sub([ops.Ny, ops.Nx],peak_inds(i));
                    df_map_i(X>col+1) = 0;
                    df_map_i(X<col-1) = 0;
                    df_map_i(Y>row+1) = 0;
                    df_map_i(Y<row-1) = 0;
                    [~, peak_inds_i] = max(df_map_i,[],"all","linear");
                    % in the rare occation that two peaks merge into one
                    if ~isempty(peak_x) && peak_x(end) == X(peak_inds_i) && ...
                            peak_y(end) == Y(peak_inds_i) && ...
                            peak_t(end) == t(k)
                        continue
                    end
                    peak_x(end+1) = X(peak_inds_i);
                    peak_y(end+1) = Y(peak_inds_i);
                    peak_inds_i = sub2ind([ops.Ny, ops.Nx], Y(peak_inds_i), X(peak_inds_i));
                    [~,peak_inds_i] = intersect(ind, peak_inds_i,'stable');
                    peak_amp(end+1) = signal_dfof(t(k),peak_inds_i);
                    peak_t(end+1) = t(k);
                end
           end
            % % Plot
            % figure
            % hold on
            % imagesc(x,y,df_map)
            % title(sprintf('ROI%03d t = %.2fs',n,ops.t(t(k))))
            % if n_peak ~= 0
            %     scatter(peak_x(end-n_peak+1:end), peak_y(end-n_peak+1:end))
            % else 
            %     scatter(peak_x(end), peak_y(end))
            % end
            % axis image

        end

        % overall
        df_map = df_map0;

        for i = 1:size(idx,1)
            x_loc = PixelList(i,1)-xLim(1)+1;
            y_loc = PixelList(i,2)-yLim(1)+1;
            df_map(y_loc, x_loc) =  max(signal_df(:,idx(i)));
        end


    end
     
    if ROI(n).Area == 1  % only one pixel
        n_cluster = n_cluster + 1 ;
        event_cluster(n_cluster).ROI = n;
        event_cluster(n_cluster).x = PixelList(:,2);
        event_cluster(n_cluster).y = PixelList(:,1);
        event_cluster(n_cluster).dff_t = ROI(n).dfft;
        event_cluster(n_cluster).dff = reshape(signal_dfof(ROI(n).dfft,idx),1,[]);
        event_cluster(n_cluster).n_spikes = length(ROI(n).dfft);
        event_cluster(n_cluster).ST = 0;

    elseif length(peak_x) == 1 % only one event
        n_cluster = n_cluster + 1 ;
        event_cluster(n_cluster).ROI = n;
        event_cluster(n_cluster).x = peak_x;
        event_cluster(n_cluster).y = peak_y;
        event_cluster(n_cluster).dff_t = peak_t;
        event_cluster(n_cluster).dff = peak_amp;
        event_cluster(n_cluster).n_spikes = 1;
        event_cluster(n_cluster).ST = 0;


    elseif ~isempty(peak_x) % more than one event and more than one pixel
        % clustering
        Z = linkage([peak_x',peak_y'],'centroid');
        T = cluster(Z, 'Cutoff', ops.cutoff, 'criterion', 'distance');
        for i = 1:max(T)
            n_cluster = n_cluster + 1 ;
            event_cluster(n_cluster).ROI = n;
            event_cluster(n_cluster).x = peak_x(T==i);
            event_cluster(n_cluster).y = peak_y(T==i);
            event_cluster(n_cluster).dff_t = peak_t(T==i);
            event_cluster(n_cluster).dff = peak_amp(T==i);
            event_cluster(n_cluster).n_spikes = nnz(T==i);
            event_cluster(n_cluster).ST = 0;

        end
    end

    if ~isempty(T)
        n_cluster_ROI = max(T);
        hold on
        imagesc(x,y,df_map)
        scatter(peak_x,peak_y,'red',"x")
        axis image
        xlim(xLim)
        ylim(yLim)
        set(gca, 'YDir','reverse')
        title(sprintf('ROI: %03d, %01d clusters', n, n_cluster_ROI))
        xlabel('X [px]')
        ylabel('Y [px]')
        % save figure
        saveas(gcf, strcat(ops.savedir_ROIpxMap,filesep, sprintf('pxMap_ROI%03d',n),'.fig'))
        saveas(gcf, strcat(ops.savedir_ROIpxMap,filesep, sprintf('pxMap_ROI%03d',n),ops.fig_format))
        clf(gcf)
    end
end
close(gcf)
% remove empty entries in structure array
event_cluster(n_cluster+1:end,:) = [];
toc

%% event statistics
stats = event_stats(event_cluster, ops);

%% save data
save_data()

%%
show_label_mask(event_cluster, ROI, ops)

%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%                                                                         %
%                            local functions                              %
%                                                                         %       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

function remove_px(k)
% keeps only pixels indexed in array k
    assignin('base','k',k);

    evalin('base','px = px(k,:);');
    evalin('base','ind = ind(k);');

    evalin('base','signal_raw = signal_raw(:,k);');

    evalin('base','signal_baseline = signal_baseline(:,k);');
    evalin('base','signal_baseline_movemean = signal_baseline_movemean(:,k);');
    
    evalin('base','signal_df = signal_df(:,k);');
    evalin('base','signal_df_movemean = signal_df_movemean(:,k);');

    evalin('base','signal_dfof = signal_dfof(:,k);');
    evalin('base','signal_dfof_movemean = signal_dfof_movemean(:,k);');
    
    % check if variable has been created
    W = evalin('base','whos');
    doesExist = ismember('signal_edge',{W(:).name});
    if doesExist
        evalin('base','signal_edge = signal_edge(:,k);');
        evalin('base','signal_edge_nve = signal_edge_nve(:,k);');
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

