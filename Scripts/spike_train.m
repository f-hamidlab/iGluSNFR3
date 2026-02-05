%% spike_train
% Performs spike train analysis for individual ROIs using MLspike algorithm
%
% DESCRIPTION:
%   Extracts and analyzes spike trains from ROI signals using the MLspike
%   (Markov chain Monte Carlo spike inference) method. Estimates spike times and
%   amplitudes from the averaged fluorescence signal within an ROI.
%
% USAGE:
%   1) spike_train(k)
%   2) spike_train([k1, k2, k3])  % multiple ROIs
%
% INPUTS:
%   - k: (int or int array) index/indices of ROI(s) to analyze
%
% OUTPUTS:
%   Updates base workspace variables:
%   - event_cluster: updated with spike train analysis results
%   - ROI: ROI data structure
%   - px: pixel data structure
%
% NOTES:
%   Accesses base workspace variables: event_cluster, ROI, px, signal_raw, ops
%   Requires MLspike toolbox and configured spike train parameters
%
% Last updated: 2026-02-03 15:30

function spike_train(k)
    % get variables from base workspace
    event_cluster = evalin('caller', 'event_cluster');
    ROI = evalin('caller', 'ROI');
    px = evalin('caller', 'px');
    signal_raw = evalin('caller', 'signal_raw');
    ops = evalin('caller', 'ops');

    ops.savedir_ROIpx = fullfile(ops.savedir,'ROI_px');
    if ~exist(ops.savedir_ROIpx, 'dir')
        mkdir(ops.savedir_ROIpx)
    end

    % remove ROI from event_cluster
    event_cluster = remove_ROI(event_cluster, k);
    
    % for each number in n
    label = vertcat(px.label);
    for i = 1:length(k)

        switch ops.ST.option
            case "MLspike"
                % get averaged signal of ROI
                idx = label == k(i);
                synapse_raw = mean(signal_raw(:,idx),2);
                
        
                % use spike train method
                % calculate new baseline
        
                synapse_baseline  = sliding_window_filter(synapse_raw, ops.baseline_percentage, ops.sl_window_ST);
                synapse_baseline = synapse_baseline - synapse_baseline(1);
                synapse_df = synapse_raw - synapse_baseline;
                normalised_signal = synapse_df ./ abs(synapse_df(1));
        
            
                % MLspike
                % [n, Ffit, ~, ~, ~, ~] = tps_mlspikes(normalised_signal,ops.par);          
                % event_cluster(end+1).dff_t = find(n)';
                % synapse_dfof = Ffit-1;
                % event_cluster(end).dff = synapse_dfof(event_cluster(end).dff_t)';
        
                ops.par.a = max(normalised_signal)-1;
                [n, Ffit, ~, ~, ~, ~] = tps_mlspikes(normalised_signal,ops.par); 
                [pks, locs] = findpeaks (Ffit-1, 'MinPeakDistance',4, 'MinPeakHeight', (ops.par.a/2));
                event_cluster(end+1).dff_t = reshape(locs,1,[]);
                synapse_dfof = Ffit-1;
                event_cluster(end).dff = synapse_dfof(event_cluster(end).dff_t)';

                % plot signal
                figure;
                sgtitle(sprintf('ROI: %03d', k(i)))
        
                subplot(2,1,1) % df
                hold on
                plot(ops.t, normalised_signal)
                plot(ops.t, Ffit,'k')
                ylabel('Normalised signal')
        
                subplot(2,1,2) % dfof
                hold on
                plot(ops.t, synapse_dfof)
                scatter(ops.t(event_cluster(end).dff_t),event_cluster(end).dff)
                ylabel('dfof')
                xlabel('Time [s]')
                
            
            case "findpeaks"
                [pks, locs] = findpeaks (ROI(k(i)).dfof, 'MinPeakDistance',4, 'MinPeakHeight', ops.ST.MinPeakHeight);
                event_cluster(end+1).dff_t = reshape(locs,1,[]);
                event_cluster(end).dff = ROI(k(i)).dfof(event_cluster(end).dff_t)';

                % plot signal
                figure;
                title(sprintf('ROI: %03d', k(i)))
                hold on
                plot(ops.t, ROI(k(i)).dfof)
                scatter(ops.t(event_cluster(end).dff_t),event_cluster(end).dff)
                ylabel('dfof')
                xlabel('Time [s]')
        end
        
        % save figure
        saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_ST',k(i)),'.fig'))
        saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_ST',k(i)),ops.fig_format))
        close(gcf)

        event_cluster(end).ROI = k(i);
        event_cluster(end).n_spikes = length(event_cluster(end).dff_t);
        event_cluster(end).x = ROI(k(i)).Centroid(1);
        event_cluster(end).y = ROI(k(i)).Centroid(2);
        event_cluster(end).x_weighted = event_cluster(end).x;
        event_cluster(end).y_weighted = event_cluster(end).y;
        event_cluster(end).ST = 1;
        event_cluster(end).dfof = ROI(k(i)).dfof;


    end

    assignin("caller","event_cluster",event_cluster)
end
