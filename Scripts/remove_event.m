%% remove_event
% Removes manually-specified event(s) from an ROI
%
% DESCRIPTION:
%   Deletes event time point(s) from an ROI's event list and re-clusters the
%   remaining events. Useful for removing detected false positives.
%
% USAGE:
%   1) remove_event(n, t)
%   2) remove_event(n, [t1, t2, t3])  % multiple times
%
% INPUTS:
%   - n: (int) index of ROI
%   - t: (numeric) time point(s) of event(s) to remove, unit: [s]
%
% OUTPUTS:
%   Updates base workspace variables:
%   - event_cluster: updated with re-clustered events for ROI n
%   - ROI: updated ROI structure with event times removed
%
% NOTES:
%   Accesses base workspace variables: event_cluster, ROI, ops, ind, signal_df, signal_dfof
%   Automatically re-clusters remaining events after removal
%
% Last updated: 2024-02-20 15:28
%   Modified to use evalin for base workspace access

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