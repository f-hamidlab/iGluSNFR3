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
    par.dt = 0.0100;
    par.F0 = [0.9 1.05];
    par.a = 0.75;  
    par.tau = 0.0629;
    par.drift.parameter = 0.02;
    par.algo.cmax = 1; 
    par.algo.nc = 50;  
    par.algo.nb = 50; 
    
    ops.par = par;
end
