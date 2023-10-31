%% spike_train_par
% defines parameter for spike train analysis
%
% USAGE: 
% 1) ops = spike_train_par(ops)
%
% INPUTS:
%     - ops: (struct) options and parameters

function ops = spike_train_par(ops)
    % parameters for spike train 
    dt = 1/ops.fs;
    par = tps_mlspikes('par');
    % - par     parameter structures - fields are:
    %           .dt     frame acquisition time
    %           .F0     baseline fluorescence (use [] to estimate it)
    %           .a      amplitude of 1 spike (in Delta F)
    %           .tau    decay time
    %           .spikerate
    %           .drift  subfield 'method' is 'state' or 'basis functions'
    par.dt = dt;
    par.F0 = [];
    par.a = 300; % 300
    t_half = 0.1; % 0.1
    par.tau = t_half/log(2);
    
    ops.par = par;
end
