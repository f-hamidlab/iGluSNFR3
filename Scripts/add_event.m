%% add_event
% Adds manually-specified event(s) to an ROI
%
% DESCRIPTION:
%   Inserts event time point(s) into an ROI's event list and re-clusters all
%   events. Useful for adding missed true positives.
%
% USAGE:
%   1) add_event(n, t)
%   2) add_event(n, [t1, t2, t3])  % multiple times
%
% INPUTS:
%   - n: (int) index of ROI
%   - t: (numeric) time point(s) of event(s) to add, unit: [s]
%
% OUTPUTS:
%   Updates base workspace variables:
%   - event_cluster: updated with re-clustered events for ROI n
%   - ROI: updated ROI structure with event times added
%
% NOTES:
%   Accesses base workspace variables: event_cluster, ROI, ops, ind, signal_df, signal_dfof
%   Automatically re-clusters all events after addition
%
% Last updated: 2024-02-20 15:28
%   Modified to use evalin for base workspace access

function add_event(n, t)
    
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

    % add time points to ROI(n).t
    ROI(n).t = unique(sort(horzcat(ROI(n).t, t')));

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
    % saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),'.fig'))
    % saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),ops.fig_format))
    % close(gcf)
   
end