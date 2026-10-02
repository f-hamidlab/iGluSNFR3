%% matching_clusters
% Matches synaptic clusters across multiple experimental conditions/trials
%
% DESCRIPTION:
%   Loads processed data from multiple experiments and matches detected events
%   based on spatial proximity (distance threshold). Clusters spatially nearby
%   events across trials using hierarchical clustering. Generates visualizations
%   of matched cluster distributions.
%
% USAGE:
%   1) matching_clusters(foldername)
%
% INPUTS:
%   - foldername: (string) root folder containing multiple 'processed_data.mat' files
%       (searches recursively through subdirectories)
%
% OUTPUTS:
%   - Figures showing:
%       .ClusterMatchingFig2_Distribution: scatter plot of event coordinates colored by cluster
%       .ClusterMatchingFig3_Histogram: histogram of cluster sizes
%       .ClusterMatchingFig4_*: additional matching statistics
%   - Figures saved to foldername with .fig and image formats
%
% NOTES:
%   Distance threshold: 3 pixels for cluster matching
%   Groups data by trial number assigned from processed_data files
%
% Last updated: 2026-02-03 15:30

function matching_clusters(foldername)
    filelist = dir(strcat(foldername,filesep,'**',filesep,'processed_data.mat'));
    N_trial = length(filelist);
    % N_trial = length(filelist)-1; % ignore the last dataset which was a different experiment setting
    edges = (1:N_trial+1); % for histcount
    delta_F_over_F = cell(1,1, N_trial);
    segmented_synapse = cell(1,1, N_trial);
    for t = 1:N_trial % ignore the last dataset which was a different experiment setting
        filename = fullfile(filelist(t).folder, filelist(t).name);
        load(filename,"event_cluster", "ops", "max_dff", "labelMask")
        [event_cluster.trial] = deal(t);
        delta_F_over_F{t} = [max_dff];
        segmented_synapse{t} = [labelMask];
        if t == 1
            event_cluster_tmp = event_cluster;
            load(filename, "mask", "first_frame")
        else
            event_cluster_tmp = [event_cluster_tmp; event_cluster];
        end
    end
    
    event_cluster = event_cluster_tmp;
    clear event_cluster_tmp
    
    % clustering
    dist_threshold = 3;  % [pixels]
    x = [event_cluster.x_weighted];
    y = [event_cluster.y_weighted];
    Z = linkage([x', y'],'centroid');
    T = cluster(Z, 'Cutoff', dist_threshold, 'criterion', 'distance');
    
    fig_handle = figure;
    gscatter(x,y,T)
    axis equal
    xlim([0 ops.Nx])
    ylim([0 ops.Ny])
    xlabel('X [px]')
    ylabel('Y [px]')
    set(gca, "YDir", "reverse")
    % save figure
    fig_name = 'ClusterMatchingFig1_Distribution';
    save_figure(fig_handle, fig_name, foldername, ops.fig_format, ops.close_fig);
    
    fig_handle = figure;
    histogram(T);
    xlabel('Group')
    ylabel('No. of event_cluster')
    % save figure
    fig_name = 'ClusterMatchingFig2_Histogram';
    save_figure(fig_handle, fig_name, foldername, ops.fig_format, ops.close_fig);

    % group clusters
    N_cluster = max(T);
    event_cluster_overall = struct([]);
    for c = 1:N_cluster
        event_cluster_shortlisted = event_cluster(T==c);
        
        total_count = 0;
        for t = 1:N_trial
            idx = [event_cluster_shortlisted.trial] == t;
            event_cluster_cell = struct2cell(event_cluster_shortlisted(idx));
            event_cluster_fields = fieldnames(event_cluster_shortlisted);
            idx = ismember(event_cluster_fields,'stim_response');
            stim_response_all = event_cluster_cell(idx,:);
            stim_response_all = cellfun(@(m) m(:)',stim_response_all,'UniformOutput',0);
            stim_response_all = horzcat(stim_response_all{:});

            stim_response_unique = unique(stim_response_all);

            event_cluster_overall(c).(sprintf("trial_%d_stim_response",t)) = stim_response_unique;
            event_cluster_overall(c).(sprintf("trial_%d_stim_response_count",t)) = length(stim_response_unique);

            total_count = total_count + event_cluster_overall(c).(sprintf("trial_%d_stim_response_count",t));

        end

        event_cluster_overall(c).stim_response_count = total_count;
        event_cluster_overall(c).stim_response_pc = total_count/(ops.n_stim * N_trial);
        event_cluster_overall(c).event_cluster_idx = find(T==c);
        event_cluster_overall(c).x_weighted = mean([event_cluster(event_cluster_overall(c).event_cluster_idx).x_weighted]);
        event_cluster_overall(c).y_weighted = mean([event_cluster(event_cluster_overall(c).event_cluster_idx).y_weighted]);

    end

    % create a figure
    fig_handle = figure;

    if ops.use_binary_mask
        [row,col] = find(mask.BW);
        buffer = 10; % pixels
        min_x = max(min(col) - buffer, 1);
        max_x = min(max(col) + buffer, ops.Nx);
        min_y = max(min(row) - buffer, 1);
        max_y = min(max(row) + buffer, ops.Ny);
    else
        min_x = 1;
        max_x = ops.Nx;
        min_y = 1;
        max_y = ops.Ny;
    end

    % plot first frame of the first trial as grayscale background, overlaid with mask.BW in red transparent foreground
    ax = subplot(2,2,1);
    imagesc(first_frame);
    colormap(ax, gray)
    if ops.use_binary_mask
        hold on
        rgbMask = cat(3, ones(size(first_frame)), zeros(size(first_frame)), zeros(size(first_frame)));
        h = image(rgbMask);
        h.AlphaData = 0.3 * mask.BW;
        title('Traced Dendrite')
    else
        title('Dendrite (no mask)')
    end
    xlabel('X [px]')
    ylabel('Y [px]')
    axis image
    xlim([min_x max_x])
    ylim([min_y max_y])
    
    % plot max of delta F over F for all trials
    ax = subplot(2,2,2);
    max_dfof = max(cell2mat(delta_F_over_F),[],3);
    max_dfof = max_dfof.*mask.BW;
    imagesc(max_dfof);
    title('Max. delta F over F')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap(ax, fire(256))
    axis image
    xlim([min_x max_x])
    ylim([min_y max_y])

    % plot segmented synapses
    ax = subplot(2,2,3);
    % plot mask only, ignore labels
    imagesc(max(cell2mat(segmented_synapse),[],3)>0);
    title('Segmented Synapses')
    xlabel('X [px]')
    ylabel('Y [px]')
    colormap(ax, gray)
    axis image
    xlim([min_x max_x])
    ylim([min_y max_y])

    % plot release probability
    ax = subplot(2,2,4);
    imagesc(first_frame);
    colormap(ax, gray)
    hold on
    % set colours for scatter points based on stim_response_pc according to the parula colormap
    cmap = parula(256);
    stim_response_pc = [event_cluster_overall.stim_response_pc];
    stim_response_pc_col = cmap(round(stim_response_pc * 255) + 1, :);
    scatter([event_cluster_overall.x_weighted], [event_cluster_overall.y_weighted], 5, stim_response_pc_col, 'filled')
    title('Release Probability')
    xlabel('X [px]')
    ylabel('Y [px]')
    axis image
    xlim([min_x max_x])
    ylim([min_y max_y])

    % save figure
    fig_name = 'ClusterMatchingFig3_Overall';
    save_figure(fig_handle, fig_name, foldername, ops.fig_format, ops.close_fig);

    filename = strcat(foldername, filesep, 'results.mat');
    save(filename,"event_cluster_overall","event_cluster");

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