%% remove_event
% remove event for ROI n at time t
%
% Last updated: 2024-02-20 15:28
%               (evalin)
%
% USAGE:
% 1) remove_event(n, t)
%
% INPUTS:
%   - n: (int) index of ROI
%   - t: (double) timepoint of event to be removed, single number or
%     array, unit: [s]

function remove_event(n, t)
    
    % get variables from base workspace
    event_cluster = evalin('caller', 'event_cluster');
    ROI = evalin('caller', 'ROI');
    ops = evalin('caller', 'ops');
    ind = evalin('caller', 'ind');
    signal_df = evalin('caller', 'signal_df');
    signal_dfof = evalin('caller', 'signal_dfof');

    % remove ROI from event_cluster
    event_cluster = remove_ROI(event_cluster, n);

    % convert time t to frame number
    t = find(ismember(ops.t,t));

    % remove time points to ROI(n).t
    ismem = ismember(ROI(n).t,t);
    ROI(n).t = ROI(n).t(~ismem);

    %% clustering
    figure
    event_cluster = cluster_events(ROI,event_cluster,n);
    close(gcf)
    assignin("caller","event_cluster",event_cluster)
    assignin("caller","ROI",ROI)

    %% plot signal for ROI
    figure;
    set(gcf,'Position',[300 300 1500 400])
    title(sprintf('ROI: %03d, Area: %d',n, ROI(n).Area))
    hold on
    plot(ops.t, ROI(n).df)
    scatter(ops.t(ROI(n).t), ROI(n).df(ROI(n).t))
    ylabel('df')
    xlabel('Time [s]')
    % save figure
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),'.fig'))
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),ops.fig_format))
    close(gcf)
   
end