%% plot_signal_pixel
% Plot intensity of a pixel over time, with subplots of 
% 1) raw signal and baseline
% 2) delta F
% 3) delta F over F, with detected peaks
% 4) rising and falling slopes, with detected peaks
%
% Pixel can be defined by one of the followings:
% 1) pixel number in structure px
%    e.g. k = 393;
% 2) linear index of pixel location
%    e.g. k = find(ind == 66817);
% 3) xy index of pixel location
%    e.g. k = find(ind==sub2ind([ops.Ny, ops.Nx], y, x));
%         where x, y are the coordinates of the pixel
%
% k can also be an array
% e.g. k = [find(ind == 6098), ...
%          find(ind==sub2ind([ops.Ny, ops.Nx],139, 337)), ... % signal
%          find(ind == 179337), ...% noise
%          find(ind == 247676)]; % spike train

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
