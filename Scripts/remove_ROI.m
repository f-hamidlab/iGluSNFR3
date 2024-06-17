%% remove_ROI
% USAGE: 
% 1) event_cluster = remove_ROI(event_cluster, k)
%
% Last updated: 2023-10-30 15:16
%               (corrected error for idx)
%
% INPUTS:
%     - event_cluster: (struct) synapse data as structure array
%     - k: (int/ int array) index of ROI

function event_cluster = remove_ROI(event_cluster, k)

    ROI_list = vertcat(event_cluster.ROI);
    idx = ismember(ROI_list, k);
    event_cluster(idx) = [];

end