%% plot_signal_pixel
% Visualizes temporal signal dynamics for individual or multiple pixels
%
% DESCRIPTION:
%   Generates multi-panel plots for selected pixels showing: raw signal with baseline,
%   delta F, delta F/F with detected peaks, and rising/falling slopes. Pixels can be
%   selected by index, linear index, or xy coordinates.
%
% USAGE:
%   1) plot_signal_pixel()  % visualizes pixels defined by variable k
%
% INPUTS:
%   - Variable k (within script): pixel selection, can be:
%       * Single pixel index: k = 393
%       * Linear pixel index: k = find(ind == 66817)
%       * XY pixel coordinates: k = find(ind == sub2ind([ops.Ny, ops.Nx], y, x))
%       * Multiple pixels: k = [idx1, idx2, idx3]
%
% OUTPUTS:
%   - Figure(s) with 4 subplots per pixel showing:
%       1. Raw signal + baseline
%       2. Delta F
%       3. Delta F/F + detected peaks
%       4. Rising/falling slopes + detected peaks
%
% NOTES:
%   Accesses base workspace variables: px, ops, signal_raw, signal_baseline, signal_df, etc.
%
% Last updated: 2026-02-03 15:30

%% Defining pixel to be plotted
k = find(ind == 2247);
% k = find(ind==sub2ind([ops.Ny, ops.Nx], 370, 57));

%% Plotting
n_plot = 4;
for i = k
    % raw and baseline
    figure;
    subplot(n_plot, 1, 1)
    title(sprintf('PxIdx: %d; SNR = %.2f; STD = %.3f',px(i).idx, px(i).SNR, px(i).STD))
    hold on
    plot(ops.t, signal_raw(:,i))
    plot(ops.t, signal_baseline(:,i))
    ylabel('Raw F')
    
    % delta F
    subplot(n_plot, 1, 2)
    plot(ops.t, signal_df(:,i))
    ylabel('df')

    % delta F over F, with detected peaks
    subplot(n_plot, 1, 3)
    hold on
    plot(ops.t, signal_dfof(:,i))
    scatter(ops.t(px(i).dff_t),px(i).dff)
    scatter(ops.t(px(i).ipt),signal_dfof(px(i).ipt,i),'rx')
    plot(ops.t, signal_dfof_movemean(:,i),'k')
    yline(3*px(i).STD)
    yline(0,'k')
    ylabel('dfof')

    % rising and falling slopes, with detected peaks
    subplot(n_plot, 1, 4)
    hold on
    plot(ops.t, signal_edge(:,i))
    plot(ops.t, signal_edge_nve(:,i))
    scatter(ops.t(px(i).slope_t),px(i).slope_pk)
    ylabel('Slope')

    xlabel('Time [s]')

end
