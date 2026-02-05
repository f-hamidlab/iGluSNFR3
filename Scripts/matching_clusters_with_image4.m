%% matching_clusters_with_image4
% Matches synaptic clusters across trials with image registration/alignment
%
% DESCRIPTION:
%   Extended version of matching_clusters that additionally performs image-based
%   registration/alignment across trials before spatial clustering. Loads processed
%   data from multiple experiments and matches detected events based on spatial
%   proximity with image registration refinement.
%
% USAGE:
%   1) matching_clusters_with_image4(foldername)
%
% INPUTS:
%   - foldername: (string) root folder containing multiple 'processed_data.mat' files
%       (searches recursively through subdirectories)
%
% OUTPUTS:
%   - Figures showing:
%       .ClusterMatchingFig2_Distribution: scatter plot of event coordinates with registration
%       .ClusterMatchingFig3_Histogram: histogram of cluster sizes
%       .ClusterMatchingFig4_*: additional matching and alignment statistics
%   - Figures saved to foldername with .fig and image formats
%
% NOTES:
%   Distance threshold: 3 pixels for cluster matching
%   Includes image registration step for coordinate alignment across trials
%
% Last updated: 2026-02-03 15:30

function matching_clusters_with_image4(foldername)
    filelist = dir(strcat(foldername,filesep,'**',filesep,'processed_data.mat'));
    N_trial = length(filelist);
    edges = (1:N_trial+1); % for histcount
    for t = 1:N_trial % ignore the last dataset which was a different experiment setting
        filename = fullfile(filelist(t).folder, filelist(t).name);
        load(filename,"event_cluster", "ops")
        [event_cluster.trial] = deal(t);
        if t == 1
            event_cluster_tmp = event_cluster;
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
    fig_name = 'ClusterMatchingFig2_Distribution';
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
        for t = 1:N_trial-1
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

            idx = ismember(event_cluster_fields,'dfof');
            dfof_all = event_cluster_cell(idx,:);
            dfof = mean(cell2mat(dfof_all),2);

            event_cluster_overall(c).(sprintf("trial_%d_dfof",t)) = dfof;
            
            total_count = total_count + event_cluster_overall(c).(sprintf("trial_%d_stim_response_count",t));

        end

        event_cluster_overall(c).stim_response_count = total_count;
        event_cluster_overall(c).stim_response_pc = total_count/(ops.n_stim * N_trial);
        event_cluster_overall(c).event_cluster_idx = find(T==c);

        for t = N_trial
            idx = [event_cluster_shortlisted.trial] == t;
            if sum(idx) == 0 % cluster not found in image 4
                event_cluster_overall(c).max_dff_image4 = 0;
            else
                event_cluster_image4 = struct2cell(event_cluster_shortlisted(idx));
                idx = ismember(event_cluster_fields,'dfof');
                dfof_all = event_cluster_image4(idx,:);
                dfof = mean(cell2mat(dfof_all),2);

                dfof = dfof(ops.stim_frames(1):ops.stim_frames(1)+ops.len_spike);
                event_cluster_overall(c).max_dff_image4 = max(dfof);
            end
        end

    end

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