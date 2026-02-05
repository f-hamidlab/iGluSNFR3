%% remove_ROI
% Removes specified ROI(s) from the event_cluster
%
% DESCRIPTION:
%   Deletes one or more ROIs from the event_cluster structure. Useful for
%   removing spurious events or artifacts detected during analysis.
%
% USAGE:
%   1) event_cluster = remove_ROI(event_cluster, k)
%   2) event_cluster = remove_ROI(event_cluster, [k1, k2, k3])  % multiple ROIs
%
% INPUTS:
%   - event_cluster: (struct) synapse data as structure array
%   - k: (int or int array) index/indices of ROI(s) to remove
%
% OUTPUTS:
%   - event_cluster: (struct) updated event_cluster with specified ROIs removed
%
% Last updated: 2023-10-30 15:16
%   Corrected error handling for idx assignment

function event_cluster = remove_ROI(event_cluster, k)

    ROI_list = vertcat(event_cluster.ROI);
    idx = ismember(ROI_list, k);
    event_cluster(idx) = [];

end