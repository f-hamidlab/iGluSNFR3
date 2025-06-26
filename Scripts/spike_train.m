%% spike_train
% use spike train analysis for ROI k
%
% USAGE: 
% 1) spike_train(k)
%
% INPUTS:
%     - k: (int/ int array) index of ROI

function spike_train(k)
    % get variables from base workspace
    event_cluster = evalin('caller', 'event_cluster');
    ROI = evalin('caller', 'ROI');
    px = evalin('caller', 'px');
    signal_raw = evalin('caller', 'signal_raw');
    ops = evalin('caller', 'ops');


    % remove ROI from event_cluster
    event_cluster = remove_ROI(event_cluster, k);
    
    % for each number in n
    label = vertcat(px.label);
    for i = 1:length(k)
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
        % save figure
        saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_ST',k(i)),'.fig'))
        saveas(gcf, strcat(ops.savedir_ROIpx,filesep, sprintf('px_ROI%03d_ST',k(i)),ops.fig_format))
        close(gcf)

        event_cluster(end).ROI = k(i);
        event_cluster(end).n_spikes = length(event_cluster(end).dff_t);
        event_cluster(end).x = ROI(k(i)).Centroid(1);
        event_cluster(end).y = ROI(k(i)).Centroid(2);
        event_cluster(end).ST = 1;


    end

    assignin("caller","event_cluster",event_cluster)
end
