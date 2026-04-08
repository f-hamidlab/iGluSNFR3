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

% List of variables to potentially save
vars_to_save = {"px", "mask", "ind", "ops", "signal_raw", "signal_df", "signal_dfof", ...
                "signal_dfof_movemean", "signal_baseline", "signal_edge", "event_cluster", "ROI", "stats"};

% Check which variables exist in caller workspace
existing_vars = {};
for i = 1:length(vars_to_save)
    var_name = vars_to_save{i};
    if evalin("caller", sprintf("exist('%s', 'var')", var_name))
        existing_vars{end+1} = var_name;
    else
        evalin("caller", sprintf("disp('  Variable ''%s'' not found - will not be saved')", var_name));
    end
end

% Build and execute save command with existing variables only
if ~isempty(existing_vars)
    % Create comma-separated list of variable names as strings
    var_strings = cellfun(@(x) sprintf('"%s"', x), existing_vars, 'UniformOutput', false);
    var_list = strjoin(var_strings, ',');
    
    % Execute save command
    disp(['Saving ' num2str(length(existing_vars)) ' variables...'])
    save_cmd = sprintf('filename = strcat(ops.savedir, filesep, ''processed_data.mat''); save(filename,%s)', var_list);
    evalin("caller", save_cmd);
    disp(['Successfully saved ' num2str(length(existing_vars)) ' variables.'])
else
    warning('No variables found to save. Creating empty processed_data.mat.')
    evalin("caller", "filename = strcat(ops.savedir, filesep, 'processed_data.mat'); save(filename, 'ops');")
end

toc

end