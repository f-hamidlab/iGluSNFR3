%% event_stats
% Calculates comprehensive statistics for detected synaptic events
%
% DESCRIPTION:
%   Computes event statistics including total event counts, amplitude distributions,
%   temporal synchronization patterns, and stimulus response metrics for evoked
%   activity experiments.
%
% USAGE:
%   1) [stats, event_cluster] = event_stats(event_cluster, ops)
%
% INPUTS:
%   - event_cluster: (struct) synapse data as structure array containing detected events
%   - ops: (struct) options and parameters including:
%       .Nt: total number of time points
%       .t: time vector
%       .savedir: directory for saving figures
%       .fig_format: figure file format
%       .experiment_type: 'spontaneous' or 'evoked'
%       .n_stim: (evoked) number of stimulation trials
%       .stim_pk_search_range: (evoked) time indices for stimulus response window
%
% OUTPUTS:
%   - stats: (struct) computed statistics including:
%       .N_syn: number of synapses/ROIs
%       .n_spikes_total: total number of events
%       .dff_all: all spike amplitudes
%       .dfft_all: temporal indices of spikes
%       .stim_response_hiscount: (evoked) response counts per stimulus
%   - event_cluster: (struct) updated event_cluster with stimulus response data (evoked)
%
% Last updated: 2026-02-03 15:30

function [stats, event_cluster] = event_stats(event_cluster, ops)

    tic;
    disp('Calculating statistics...')
    stats.N_syn = size(event_cluster,1);
    % number of synapse
    fprintf('Number of synapse = %d\n', stats.N_syn);
    
    % total number of events
    stats.n_spikes_total = sum(vertcat(event_cluster.n_spikes));
    fprintf('Total number of events = %d\n', stats.n_spikes_total);
    
    % intensity of events
    event_cluster_cell = struct2cell(event_cluster);
    event_cluster_fields = fieldnames(event_cluster);
    idx = ismember(event_cluster_fields,'dff');
    dff_all = event_cluster_cell(idx,:);
    dff_all = cellfun(@(m)m(1,:),dff_all,'UniformOutput',0);
    stats.dff_all = horzcat(dff_all{:});
    fprintf('Spike amplitude (dfof) = %.3f \x00B1 %.3f \n', mean(stats.dff_all), std(stats.dff_all));
    
    % synapse firing together
    idx = ismember(event_cluster_fields,'dff_t');
    dfft_all = event_cluster_cell(idx,:);
    dfft_all = cellfun(@(m)m(1,:),dfft_all,'UniformOutput',0);
    stats.dfft_all = horzcat(dfft_all{:});
    edges = (0:ops.Nt);
    [N,~] = histcounts(stats.dfft_all,edges);
    fig_handle = figure;
    plot(ops.t, N)
    ylabel('Count')
    xlabel('Time [s]')
    title('Number of synapses firing together')
    % save figure
    saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses.fig'))
    saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses',ops.fig_format))
    close(fig_handle)

    % moving sum
    kernal = ones(5,1);
    N_w = conv(N, kernal, 'same');
    fig_handle = figure;
    plot(ops.t, N_w)
    ylabel('Count')
    xlabel('Time [s]')
    title('Number of synapses firing together (moving sum)')
    % save figure
    saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses_sum.fig'))
    saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses_sum',ops.fig_format))
    close(fig_handle)
    toc

    if ops.experiment_type == "evoked"
        % get ratio of evoked response to total number of stimulation
        for n = 1:size(event_cluster, 1)
        
            [~,I] = max(event_cluster(n).dfof(ops.stim_pk_search_range),[],2);
            
            % get different column from each row
            I = bsxfun(@eq, cumsum(ones(size(ops.stim_pk_search_range)), 2), I);
            event_cluster(n).stim_max_t = sum(ops.stim_pk_search_range.*I, 2);
    
            event_cluster(n).stim_max_dff = event_cluster(n).dfof(event_cluster(n).stim_max_t);
    
            % intersect is not an accurate measure
            % event_cluster(n).t_intersect = intersect(event_cluster(n).dff_t, event_cluster(n).stim_max_t);
                
            [~,Locb] = ismember(event_cluster(n).dff_t, ops.stim_pk_search_range);
            [row,~] = ind2sub(size(ops.stim_pk_search_range), Locb);
              
            row = unique(row);
            row = row(row~=0);
            event_cluster(n).stim_response = row;
            event_cluster(n).stim_response_count = length(row);
            event_cluster(n).stim_response_pc = length(row)/ops.n_stim;
                
        end

        % percentage of event cluster responding to each stim
        event_cluster_cell = struct2cell(event_cluster);
        event_cluster_fields = fieldnames(event_cluster);
        idx = ismember(event_cluster_fields,'stim_response');
        stim_response_all = event_cluster_cell(idx,:);
        stim_response_all = cellfun(@(m) m(:)',stim_response_all,'UniformOutput',0);
        stim_response_all = horzcat(stim_response_all{:});
        edges = (1:ops.n_stim+1);
        [stats.stim_response_hiscount,~] = histcounts(stim_response_all,edges);
        fig_handle = figure;
        plot(stats.stim_response_hiscount)
        ylabel('Count')
        xlabel('Stimulation')
        title('Number of synapses firing together')
        % save figure
        saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig9_StimulationResponse.fig'))
        saveas(fig_handle, strcat(ops.savedir,filesep, 'Fig9_StimulationResponse',ops.fig_format))
        close(fig_handle)
    end

end