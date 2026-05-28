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
    filelist = dir(strcat(foldername, filesep, '**', filesep, 'processed_data.mat'));
    if isempty(filelist)
        warning('No processed_data.mat files found in %s', foldername);
        return
    end

    [group_size, group_records] = build_group_records(filelist);
    if isempty(group_records)
        warning('No valid processed_data.mat files found for matching in %s', foldername);
        return
    end

    group_keys = unique({group_records.group_key}, 'stable');
    for g = 1:numel(group_keys)
        group_key = group_keys{g};
        idx_group = strcmp({group_records.group_key}, group_key);
        group_data = group_records(idx_group);
        [~, sort_idx] = sort([group_data.rec_num]);
        group_data = group_data(sort_idx);

        event_cluster_tmp = struct([]);
        valid_trial_count = 0;
        for t = 1:numel(group_data)
            s = load(group_data(t).filepath, 'event_cluster', 'ops');
            if ~isfield(s, 'event_cluster') || isempty(s.event_cluster)
                disp(['No event_cluster variable found in ', group_data(t).filepath, ' - skipping'])
                continue
            end

            valid_trial_count = valid_trial_count + 1;
            [s.event_cluster.trial] = deal(valid_trial_count);
            if isempty(event_cluster_tmp)
                event_cluster_tmp = s.event_cluster;
            else
                event_cluster_tmp = [event_cluster_tmp; s.event_cluster];
            end
            ops = s.ops;
        end

        if isempty(event_cluster_tmp)
            warning('No non-empty event_cluster found for group %s in %s', group_key, foldername);
            continue
        end

        event_cluster = event_cluster_tmp;
        N_trial = valid_trial_count;

        % clustering
        dist_threshold = 3;  % [pixels]
        x = [event_cluster.x_weighted];
        y = [event_cluster.y_weighted];
        Z = linkage([x', y'], 'centroid');
        T = cluster(Z, 'Cutoff', dist_threshold, 'criterion', 'distance');

        fig_suffix = '';
        out_suffix = '';
        if numel(group_keys) > 1
            fig_suffix = ['_', sanitize_label(group_key)];
            out_suffix = ['_', sanitize_label(group_key)];
        end

        fig_handle = figure;
        gscatter(x, y, T)
        axis equal
        xlim([0 ops.Nx])
        ylim([0 ops.Ny])
        xlabel('X [px]')
        ylabel('Y [px]')
        set(gca, "YDir", "reverse")
        fig_name = ['ClusterMatchingFig2_Distribution', fig_suffix];
        save_figure(fig_handle, fig_name, foldername, ops.fig_format, ops.close_fig);

        fig_handle = figure;
        histogram(T);
        xlabel('Group')
        ylabel('No. of event_cluster')
        fig_name = ['ClusterMatchingFig2_Histogram', fig_suffix];
        save_figure(fig_handle, fig_name, foldername, ops.fig_format, ops.close_fig);

        N_cluster = max(T);
        event_cluster_overall = struct([]);
        for c = 1:N_cluster
            event_cluster_shortlisted = event_cluster(T == c);

            total_count = 0;
            for t = 1:N_trial-1
                idx = [event_cluster_shortlisted.trial] == t;
                event_cluster_cell = struct2cell(event_cluster_shortlisted(idx));
                event_cluster_fields = fieldnames(event_cluster_shortlisted);
                idx = ismember(event_cluster_fields, 'stim_response');
                stim_response_all = event_cluster_cell(idx, :);
                stim_response_all = cellfun(@(m) m(:)', stim_response_all, 'UniformOutput', 0);
                stim_response_all = horzcat(stim_response_all{:});

                stim_response_unique = unique(stim_response_all);

                event_cluster_overall(c).(sprintf("trial_%d_stim_response", t)) = stim_response_unique;
                event_cluster_overall(c).(sprintf("trial_%d_stim_response_count", t)) = length(stim_response_unique);

                idx = ismember(event_cluster_fields, 'dfof');
                dfof_all = event_cluster_cell(idx, :);
                dfof = mean(cell2mat(dfof_all), 2);

                event_cluster_overall(c).(sprintf("trial_%d_dfof", t)) = dfof;

                total_count = total_count + event_cluster_overall(c).(sprintf("trial_%d_stim_response_count", t));
            end

            event_cluster_overall(c).stim_response_count = total_count;
            event_cluster_overall(c).stim_response_pc = total_count / (ops.n_stim * N_trial);
            event_cluster_overall(c).event_cluster_idx = find(T == c);
            event_cluster_overall(c).matching_group_size = group_size;
            event_cluster_overall(c).matching_group_key = group_key;

            idx = [event_cluster_shortlisted.trial] == N_trial;
            if sum(idx) == 0
                event_cluster_overall(c).max_dff_image4 = 0;
            else
                event_cluster_fields = fieldnames(event_cluster_shortlisted);
                event_cluster_image4 = struct2cell(event_cluster_shortlisted(idx));
                idx = ismember(event_cluster_fields, 'dfof');
                dfof_all = event_cluster_image4(idx, :);
                dfof = mean(cell2mat(dfof_all), 2);

                dfof = dfof(ops.stim_frames(1):ops.stim_frames(1)+ops.len_spike);
                event_cluster_overall(c).max_dff_image4 = max(dfof);
            end
        end

        if numel(group_keys) == 1
            filename = fullfile(foldername, 'results.mat');
        else
            filename = fullfile(foldername, ['results_group', out_suffix, '.mat']);
        end
        save(filename, 'event_cluster_overall', 'event_cluster');
    end

end

function [group_size, records] = build_group_records(filelist)
    records = struct('filepath', {}, 'cell_token', {}, 'rec_num', {}, 'group_start', {}, 'group_key', {});
    group_size = [];

    for i = 1:numel(filelist)
        filepath = fullfile(filelist(i).folder, filelist(i).name);
        [cell_token, rec_num] = parse_cell_and_rec(filepath);
        if isempty(cell_token) || isnan(rec_num)
            continue
        end

        s = load(filepath, 'ops');
        if isempty(group_size) && isfield(s, 'ops') && isfield(s.ops, 'binary_mask_group_size')
            group_size = s.ops.binary_mask_group_size;
        end

        records(end+1).filepath = filepath; %#ok<AGROW>
        records(end).cell_token = cell_token;
        records(end).rec_num = rec_num;
    end

    if isempty(group_size) || ~isscalar(group_size) || ~isfinite(group_size) || group_size < 1
        group_size = numel(records);
    end
    group_size = max(1, round(group_size));

    for i = 1:numel(records)
        group_start = floor((records(i).rec_num - 1) / group_size) * group_size + 1;
        records(i).group_start = group_start;
        records(i).group_key = sprintf('%s_%d', records(i).cell_token, group_start);
    end
end

function [cell_token, rec_num] = parse_cell_and_rec(path_str)
    cell_token = '';
    rec_num = NaN;

    token = regexp(path_str, '(Cell\d+)_(\d+)', 'tokens');
    if isempty(token)
        return
    end
    last_token = token{end};
    cell_token = last_token{1};
    rec_num = str2double(last_token{2});
end

function out = sanitize_label(in)
    out = regexprep(in, '[^A-Za-z0-9_\-]', '_');
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