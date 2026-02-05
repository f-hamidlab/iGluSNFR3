%% rasterplot
% Generates spike raster plot for all detected events and optimizes epsilon parameter
%
% DESCRIPTION:
%   Creates a raster plot visualization of spike times across all detected events/ROIs.
%   Uses DBSCAN clustering on spike data to find the optimal epsilon parameter that
%   minimizes background noise while preserving spike clusters.
%
% USAGE:
%   1) rasterplot(ops, event_cluster)
%
% INPUTS:
%   - ops: (struct) options and parameters including:
%       .fs: sampling frequency
%       .Nt: total number of time points
%       .spontaneous.peakWidth: width of spike peak
%       .savedir: directory for saving figures
%       .fig_format: figure file format
%       .close_fig: whether to close figures after saving
%   - event_cluster: (struct) synaptic event data
%
% OUTPUTS:
%   - Figure showing noise fraction vs epsilon parameter
%   - Figures saved to ops.savedir
%
% NOTES:
%   Tests epsilon values from 1 to 10 with 0.5 increments using DBSCAN clustering
%
% Last updated: 2026-02-03 15:30

function rasterplot(ops,event_cluster)
    % getting data
    spikeTimes_fields = fieldnames(event_cluster);
    spikeTimes = struct2cell(event_cluster);
    idx = ismember(spikeTimes_fields,'dff_t');
    spikeTimes = spikeTimes(idx,:);
    spikeTimes = cellfun(@(x) x./ops.fs, spikeTimes, 'UniformOutput',0);
    
    N_cluster = length(event_cluster);
    signal = false(length(event_cluster),ops.Nt);
    for c = 1:N_cluster
        t = event_cluster(c).dff_t;
        signal(c, t) = true;
    end

    % clustering
    % moving max
    signal = movmax(signal,round(ops.spontaneous.peakWidth/2),2);
        
    minpts = 5;
    epsilon_start = 1;
    epsilon_step = 0.5;
    epsilon_end = 10;
    epsilon = epsilon_start: epsilon_step :epsilon_end;
    noise_fr = zeros(1,length(epsilon));
    for i = 1:length(epsilon)
        idx = dbscan(single(signal),epsilon(i),minpts);
        noise_fr(i) = 1-(sum(idx == -1)/N_cluster);
    end

    figure_handle = figure; 
    plot(noise_fr,epsilon)
    ylabel('epsilon')
    xlabel('noise_fr')

    % save figure
    fig_name = 'dbscan';
    saveas(figure_handle, fullfile(ops.savedir, [fig_name,'.fig']))
    saveas(figure_handle, fullfile(ops.savedir, [fig_name, ops.fig_format]))
    close(figure_handle)

    [~, idx_of_result] = knee_pt(noise_fr);
    eps = epsilon(idx_of_result);
    idx = dbscan(single(signal),eps,minpts);
    if all(idx == -1)
        eps = epsilon(idx_of_result+1);
        idx = dbscan(single(signal),eps,minpts);
    end
    
    [~,idx] = sort(idx);
    spikeTimes = spikeTimes(idx);

    % plotting
    figure_handle = figure;
    LineFormat.Color = 'b';
    plotSpikeRaster(spikeTimes,'SpikeDuration',ops.spontaneous.peakWidth/ops.fs,'LineFormat',LineFormat);
    ylabel('Event cluster')
    xlabel('Time [s]')
    xlim([ops.t(1) ops.t(end)])
    
    % save figure
    fig_name = 'rasterplot';
    saveas(figure_handle, fullfile(ops.savedir, [fig_name,'.fig']))
    saveas(figure_handle, fullfile(ops.savedir, [fig_name, ops.fig_format]))
    close(figure_handle)


end