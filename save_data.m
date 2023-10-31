%% save_data
% saves results as .mat for future analysis
%
% USAGE:
% 1) save_data()
% 

function save_data()

% save data
tic;
disp('Saving data...')
evalin("base","filename = strcat(ops.savedir, filesep, 'processed_data.mat');")
evalin("base","save(filename,""px"",""ind"",""ops"",""signal_raw"",""signal_df"",""signal_dfof"",""signal_dfof_movemean"",""signal_baseline"",""signal_edge"",""signal_edge_nve"",""event_cluster"", ""ROI"",""stats"")")
toc

end