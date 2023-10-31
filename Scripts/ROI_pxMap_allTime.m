%% ROI_pxMap_allTime
% plots map of delta F for ROI n at each time point of event
%
% USAGE:
% 1) ROI_pxMap_allTime(n)
%
% INPUTS:
%   - n: (int) index of ROI

function ROI_pxMap_allTime(n)

    % get variables from base workspace
    ROI = evalin('base', 'ROI');
    ops = evalin('base', 'ops');
    ind = evalin('base', 'ind');
    signal_df = evalin('base', 'signal_df');
    signal_dfof = evalin('base', 'signal_dfof');

    % for each ROI
    fprintf('ROI%03d\n',n)
    
    % find stats of corresponding pixels
    [~,idx] = intersect(ind, ROI(n).PixelIdxList,'stable');

    t = ROI(n).t;
    
    % get location of each spike
    PixelList = ROI(n).PixelList;
    xLim = [min(PixelList(:,1))-5, max(PixelList(:,1))+5];
    yLim = [min(PixelList(:,2))-5, max(PixelList(:,2))+5];
    x = xLim(1):xLim(2);
    y = yLim(1):yLim(2);
    [X,Y] = meshgrid(x,y);

    df_map0 = zeros(yLim(2)-yLim(1)+1,xLim(2)-xLim(1)+1);

    if ROI(n).Area == 1 % only one pixel
        fprintf("ROI%03d consists of one pixel only. Skipping plot...", n);
 
    else % calculate location of event
        peak_x = [];
        peak_y = [];
        peak_amp = [];
        for k = 1: length(t)
            df_map = df_map0;
            
            for i = 1:size(idx,1)
                x_loc = PixelList(i,1)-xLim(1)+1;
                y_loc = PixelList(i,2)-yLim(1)+1;
                df_map(y_loc, x_loc) = signal_df(t(k),idx(i));
            end
            df_map(df_map<0)=0;
            z = medfilt2(df_map);
            peak_inds = findpeaks_2d(X, z, ops.MinPeakHeight);
            
            n_peak = sum(peak_inds,'all');

            if n_peak <= 1 % if there is only one spatial peak
                signal_dfof1 = signal_dfof(t(k),idx(:));
                peak_amp(end+1) = max(signal_dfof1);
        
                signal_df1 = signal_df(t(k),idx(:));
                signal_df1(signal_df1<max(signal_df1)/2) = 0;
                
                % weighted average
                peak_x(end+1) = sum(PixelList(:,1).*signal_df1')/sum(signal_df1);
                peak_y(end+1) = sum(PixelList(:,2).*signal_df1')/sum(signal_df1);

            elseif n_peak ~= 0 % more than one peak
                peak_inds = sub2ind([ops.Ny, ops.Nx], Y(peak_inds), X(peak_inds));
                for i = 1:n_peak
                    df_map_i = df_map;
                    [row, col] = ind2sub([ops.Ny, ops.Nx],peak_inds(i));
                    df_map_i(X>col+1) = 0;
                    df_map_i(X<col-1) = 0;
                    df_map_i(Y>row+1) = 0;
                    df_map_i(Y<row-1) = 0;
                    [~, peak_inds_i] = max(df_map_i,[],"all","linear");
                    peak_x(end+1) = X(peak_inds_i);
                    peak_y(end+1) = Y(peak_inds_i);
                end

            end

            % Plot Figure
            figure
            hold on
            imagesc(x,y,df_map)
            title(sprintf('ROI%03d t = %.2fs',n,ops.t(t(k))))
            if n_peak ~= 0
                scatter(peak_x(end-n_peak+1:end), peak_y(end-n_peak+1:end),'red',"x")
            else 
                scatter(peak_x(end), peak_y(end),'red',"x")
            end
            axis image
            set(gca, 'YDir','reverse')
            xlabel('X [px]')
            ylabel('Y [px]')


        end
    end
end
