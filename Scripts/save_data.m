%% save_data
% Saves processed analysis results as a MATLAB .mat file
%
% DESCRIPTION:
%   Packages all analysis results including pixel data, ROIs, signals, and
%   computed statistics into a single .mat file for downstream analysis and archival.
%
% USAGE:
%   1) save_data()
%
% INPUTS:
%   None (accesses variables from caller workspace)
%
% OUTPUTS:
%   - processed_data.mat file containing:
%       .px: pixel data and analysis results
%       .mask: spatial mask information
%       .ind: pixel indices
%       .ops: analysis parameters
%       .signal_raw: raw fluorescence signals
%       .signal_df: delta F signals
%       .signal_dfof: delta F/F signals
%       .signal_dfof_movemean: moving average of delta F/F
%       .signal_baseline: baseline fluorescence
%       .signal_edge: Gaussian edge detection output
%       .event_cluster: detected synaptic events
%       .ROI: region of interest definitions
%       .stats: computed event statistics
%
% Last updated: 2026-02-03 15:30

function save_data()

% save data
tic;
disp('Saving data...')
evalin("caller","filename = strcat(ops.savedir, filesep, 'processed_data.mat');")
evalin("caller","save(filename,""px"",""mask"",""ind"",""ops"",""signal_raw"",""signal_df"",""signal_dfof"",""signal_dfof_movemean"",""signal_baseline"",""signal_edge"",""event_cluster"", ""ROI"",""stats"")")
toc

end