%% cluster_events
% cluster events in ROI
%
% Last updated: 2025-02-04 14:37
%               (fixed PixelList)
%
% USAGE:
% 1) event_cluster = cluster_events(ROI,event_cluster,n)
%
% INPUTS:
%     - ROI: (struct) ROI data as structure array
%     - event_cluster: (struct) synapse data as structure array
%     - n: (int) index of ROI
%
% OUTPUTS:
%     - event_cluster: (struct) synapse data as structure array


function event_cluster = cluster_events(ROI,event_cluster,n)
    % get variables from base workspace    
    ops = evalin('caller', 'ops');
    ind = evalin('caller', 'ind');
    signal_df = evalin('caller', 'signal_df');
    signal_dfof = evalin('caller', 'signal_dfof');

    fprintf('ROI%03d\n',n)
    
    % find stats of corresponding pixels
    [~,idx] = intersect(ind, ROI(n).PixelIdxList,'stable');

    t = ROI(n).t;
    
    % get location of each spike
    [row,col] = ind2sub([ops.Ny, ops.Nx],ind(idx));
    PixelList = [col,row];

    xLim = [min(PixelList(:,1))-5, max(PixelList(:,1))+5];
    yLim = [min(PixelList(:,2))-5, max(PixelList(:,2))+5];
    x = xLim(1):xLim(2);
    y = yLim(1):yLim(2);
    [X,Y] = meshgrid(x,y);

    df_map0 = zeros(yLim(2)-yLim(1)+1,xLim(2)-xLim(1)+1);
    T = [];

    if ROI(n).Area == 1 % only one pixel
        peak_x = PixelList(:,2);
 
    else % calculate location of event
        peak_x = [];
        peak_y = [];
        peak_amp = [];
        peak_t = [];
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

                peak_t(end+1) = t(k);

            elseif n_peak ~= 0 % more than one peak
                peak_x(end+1:end+n_peak) = X(peak_inds);
                peak_y(end+1:end+n_peak) = Y(peak_inds);
                peak_inds = sub2ind([ops.Ny, ops.Nx], Y(peak_inds), X(peak_inds));
                for i = 1:n_peak
                    [~,peak_inds_i] = intersect(ind, peak_inds(i),'stable');
                    if ~isempty(peak_inds_i)
                        peak_amp(end+1) = signal_dfof(t(k),peak_inds_i);
                        peak_t(end+1) = t(k);
                    else %
                        df_map_i = df_map;
                        [row, col] = ind2sub([ops.Ny, ops.Nx],peak_inds(i));
                        df_map_i(X>col+1) = 0;
                        df_map_i(X<col-1) = 0;
                        df_map_i(Y>row+1) = 0;
                        df_map_i(Y<row-1) = 0;
                        [~, peak_inds_i] = max(df_map_i,[],"all","linear");
                        peak_inds_i = sub2ind([ops.Ny, ops.Nx], Y(peak_inds_i), X(peak_inds_i));
                        [~,peak_inds_i] = intersect(ind, peak_inds_i,'stable');
                        peak_amp(end+1) = signal_dfof(t(k),peak_inds_i);
                        peak_t(end+1) = t(k);
                    end
                end
           end
            % % Plot
            % figure
            % hold on
            % imagesc(x,y,df_map)
            % title(sprintf('ROI%03d t = %.2fs',n,ops.t(t(k))))
            % if n_peak ~= 0
            %     scatter(peak_x(end-n_peak+1:end), peak_y(end-n_peak+1:end))
            % else 
            %     scatter(peak_x(end), peak_y(end))
            % end
            % axis image

        end

        % overall
        df_map = df_map0;

        for i = 1:size(idx,1)
            x_loc = PixelList(i,1)-xLim(1)+1;
            y_loc = PixelList(i,2)-yLim(1)+1;
            df_map(y_loc, x_loc) =  max(signal_df(:,idx(i)));
        end
    end
     
    if ROI(n).Area == 1  % only one pixel
        event_cluster(end+1).ROI = n;
        event_cluster(end).x = PixelList(:,1);
        event_cluster(end).y = PixelList(:,2);
        event_cluster(end).dff_t = ROI(n).t;
        event_cluster(end).dff = signal_dfof(ROI(n).t,idx);
        event_cluster(end).n_spikes = length(ROI(n).t);
        event_cluster(end).ST = 0;

    elseif isscalar(peak_x) % only one event
        event_cluster(end+1).ROI = n;
        event_cluster(end).x = peak_x;
        event_cluster(end).y = peak_y;
        event_cluster(end).dff_t = peak_t;
        event_cluster(end).dff = peak_amp;
        event_cluster(end).n_spikes = 1;
        event_cluster(end).ST = 0;


    elseif ~isempty(peak_x) % more than one event and more than one pixel
        % clustering
        Z = linkage([peak_x',peak_y'],'centroid');
        T = cluster(Z, 'Cutoff', ops.cutoff, 'criterion', 'distance');
        
        for i = 1:max(T)
            event_cluster(end+1).ROI = n;
            event_cluster(end).ST = 0;

            x_tmp = peak_x(T==i);
            y_tmp = peak_y(T==i);
            dff_t_tmp = peak_t(T==i);
            dff_tmp = peak_amp(T==i);

            % remove repeated spikes
            [event_cluster(end).dff_t,ia,~] = unique(dff_t_tmp,'legacy');
            event_cluster(end).n_spikes = length(ia);
            
            x_weighted = zeros(1,length(ia));
            y_weighted = zeros(1,length(ia));
            dff_weighted = zeros(1,length(ia));

            for k = 1:event_cluster(end).n_spikes
                i_k = find(dff_t_tmp == event_cluster(end).dff_t(k),1);
                
                % weighted average
                x_weighted(k) = sum(x_tmp(i_k).*dff_tmp(i_k))/sum(dff_tmp(i_k));
                y_weighted(k) = sum(y_tmp(i_k).*dff_tmp(i_k))/sum(dff_tmp(i_k));
                dff_weighted(k) = sum(dff_tmp(i_k).*dff_tmp(i_k))/sum(dff_tmp(i_k));
            end
            event_cluster(end).x = x_weighted;
            event_cluster(end).y = y_weighted;
            event_cluster(end).dff = dff_weighted;
        end
    end

    if ~isempty(T)
        n_cluster_ROI = max(T);
        hold on
        imagesc(x,y,df_map)
        scatter(peak_x,peak_y,'red',"x")
        axis image
        xlim(xLim)
        ylim(yLim)
        set(gca, 'YDir','reverse')
        title(sprintf('ROI: %03d, %01d clusters', n, n_cluster_ROI))
        xlabel('X [px]')
        ylabel('Y [px]')
        % save figure
        saveas(gcf, strcat(ops.savedir_ROIpxMap,filesep, sprintf('pxMap_ROI%03d',n),'.fig'))
        saveas(gcf, strcat(ops.savedir_ROIpxMap,filesep, sprintf('pxMap_ROI%03d',n),ops.fig_format))
        clf(gcf)
    end
end
