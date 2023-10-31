%% event_stats
% calculates event statistics
%
% USAGE: 
% 1) stats = event_stats(event_cluster, ops)
%
% INPUTS:
%     - event_cluster: (struct) synapse data as structure array
%     - ops: (struct) options and parameters

function stats = event_stats(event_cluster, ops)

    tic;
    disp('Calculating statistics...')
    stats.N_syn = size(event_cluster,1);
    % number of synapse
    fprintf('Number of synapse = %d\n', stats.N_syn);
    
    % total number of events
    stats.n_spikes_total = sum(vertcat(event_cluster.n_spikes));
    fprintf('Total number of events = %d\n', stats.n_spikes_total);
    
    % intensity of events
    dff_all = struct2cell(event_cluster);
    fields = fieldnames(event_cluster);
    idx = find(ismember(fields,'dff'));
    dff_all = dff_all(idx,:);
    dff_all = cellfun(@(m)m(1,:),dff_all,'UniformOutput',0);
    stats.dff_all = horzcat(dff_all{:});
    fprintf('Spike amplitude (dfof) = %.3f \x00B1 %.3f \n', mean(stats.dff_all), std(stats.dff_all));
    
    % synapse firing together
    dfft_all = struct2cell(event_cluster);
    fields = fieldnames(event_cluster);
    idx = find(ismember(fields,'dff_t'));
    dfft_all = dfft_all(idx,:);
    dfft_all = cellfun(@(m)m(1,:),dfft_all,'UniformOutput',0);
    stats.dfft_all = horzcat(dfft_all{:});
    edges = (0:ops.Nt);
    [N,~] = histcounts(stats.dfft_all,edges);
    figure;
    plot(ops.t, N)
    ylabel('Count')
    xlabel('Time [s]')
    title('Number of synapses firing together')
    % save figure
    saveas(gcf, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses.fig'))
    saveas(gcf, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses',ops.fig_format))
    close(gcf)
    % moving sum
    kernal = ones(5,1);
    N_w = conv(N, kernal, 'same');
    figure;
    plot(ops.t, N_w)
    ylabel('Count')
    xlabel('Time [s]')
    title('Number of synapses firing together (moving sum)')
    % save figure
    saveas(gcf, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses_sum.fig'))
    saveas(gcf, strcat(ops.savedir,filesep, 'Fig8_SynchroniousSynapses_sum',ops.fig_format))
    close(gcf)
    toc

end
