%% show_label_mask_with_text
% Visualizes ROI label mask with ROI numbers labeled on centroids
%
% DESCRIPTION:
%   Creates and displays a labeled image showing all ROIs with their centroid
%   positions. ROI numbers are overlaid on the image at each centroid location.
%   Each unique ROI is assigned a distinct color.
%
% USAGE:
%   1) show_label_mask_with_text(event_cluster, ROI, ops)
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
%   - Figures saved: Fig_LabelMask_Text.fig and Fig_LabelMask_Text.[format]
%   - No variables returned
%
% Last updated: 2026-02-03 15:30

function show_label_mask_with_text(event_cluster, ROI, ops)

    labelMask = zeros(ops.Ny, ops.Nx);
    ROI_list = vertcat(event_cluster.ROI);
    ROI_list = unique(ROI_list);

    for i = 1:length(ROI_list)
        ROI_num = ROI_list(i);
        labelMask(ROI(ROI_num).PixelIdxList) = ROI_num;
    end
    
    s = regionprops(labelMask, {'Centroid','PixelIdxList'});
    
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
    
    
hold on
for k = 1:numel(s)
    centroid = s(k).Centroid;
    text(centroid(1), centroid(2), sprintf('%d', k),'Color','r','FontSize',7);
end
hold off

    
    % save figure
    saveas(gcf, strcat(ops.savedir,filesep,'Fig_LabelMask_text.fig'))
    saveas(gcf, strcat(ops.savedir,filesep,'Fig_LabelMask_text',ops.fig_format))
    % close(gcf)

end