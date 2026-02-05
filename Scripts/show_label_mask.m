%% show_label_mask
% Visualizes ROI label mask with color-coded regions
%
% DESCRIPTION:
%   Creates and displays a labeled image showing all ROIs detected in the
%   event_cluster. Each unique ROI is assigned a distinct color. Generates
%   both .fig and image file outputs.
%
% USAGE:
%   1) show_label_mask(event_cluster, ROI, ops)
%
% INPUTS:
%   - event_cluster: (struct) synapse data as structure array
%   - ROI: (struct) ROI data as structure array
%   - ops: (struct) options and parameters including:
%       .Ny, .Nx: image dimensions
%       .savedir: directory for saving figures
%       .fig_format: figure file format
%
% OUTPUTS:
%   - Figures saved: Fig_LabelMask.fig and Fig_LabelMask.[format]
%   - No variables returned
%
% Last updated: 2026-02-03 15:30

function show_label_mask(event_cluster, ROI, ops)

    labelMask = zeros(ops.Ny, ops.Nx);
    ROI_list = vertcat(event_cluster.ROI);
    ROI_list = unique(ROI_list);

    for i = 1:length(ROI_list)
        ROI_num = ROI_list(i);
        labelMask(ROI(ROI_num).PixelIdxList) = ROI_num;
    end
    
    % Plot label mask
    figure;
    imagesc(labelMask)
    title('Label Mask')
    xlabel('X [px]')
    ylabel('Y [px]')
    
    cmap = parula;
    cmap(1,:)=[0,0,0];
    colormap(cmap)
    
    % colormap parula
    axis image
    
    % save figure
    saveas(gcf, strcat(ops.savedir,filesep,'Fig_LabelMask.fig'))
    saveas(gcf, strcat(ops.savedir,filesep,'Fig_LabelMask',ops.fig_format))
    % close(gcf)

end