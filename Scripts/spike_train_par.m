%% spike_train_par
% Configures parameters for MLspike spike train analysis
%
% DESCRIPTION:
%   Sets up parameters for the MLspike (Markov chain Monte Carlo spike inference)
%   algorithm including baseline fluorescence, spike amplitude, decay time, and
%   drift correction settings.
%
% USAGE:
%   1) ops = spike_train_par(ops)
%
% INPUTS:
%   - ops: (struct) options and parameters, must contain .fs (sampling frequency)
%
% OUTPUTS:
%   - ops: (struct) updated with .par field containing MLspike parameters:
%       .par.dt: frame acquisition time
%       .par.F0: baseline fluorescence range [min, max]
%       .par.a: spike amplitude in Delta F
%       .par.tau: exponential decay time constant
%       .par.spikerate: baseline spike rate
%       .par.drift: drift correction parameters
%       .par.algo: algorithm settings (cmax, nc, nb)
%
% Last updated: 2026-02-03 15:30

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
    par.F0 = [0.9 1.05];
    par.a = 0.3;  % 0.75
    par.tau = 0.0629; % 0.0629
    par.drift.parameter = 0.02;
    par.algo.cmax = 1; 
    par.algo.nc = 50;  
    par.algo.nb = 50; 

   
    ops.par = par;
end
