%% remove_event
% remove event for ROI n at time t
%
% Last updated: 2023-10-30 15:17
%               (corrected error for horzcat)
%
% USAGE:
% 1) remove_event(n, t)
%
% INPUTS:
%   - n: (int) index of ROI
%   - t: (double) timepoint of event to be removed, single number or
%     array, unit: [s]

function remove_event(n, t)
    
    % get variables from base workspace
    event_cluster = evalin('base', 'event_cluster');
    ROI = evalin('base', 'ROI');
    ops = evalin('base', 'ops');
    ind = evalin('base', 'ind');
    signal_df = evalin('base', 'signal_df');
    signal_dfof = evalin('base', 'signal_dfof');

    % remove ROI from event_cluster
    event_cluster = remove_ROI(event_cluster, n);

    % convert time t to frame number
    t = find(ismember(ops.t,t));

    % remove time points to ROI(n).t
    ismem = ismember(ROI(n).t,t);
    ROI(n).t = ROI(n).t(~ismem);

    %% clustering
    figure

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
        event_cluster(end).x = PixelList(:,2);
        event_cluster(end).y = PixelList(:,1);
        event_cluster(end).dff_t = ROI(n).t;
        event_cluster(end).dff = signal_dfof(ROI(n).t,idx);
        event_cluster(end).n_spikes = length(ROI(n).t);
        event_cluster(end).ST = 0;

    elseif length(peak_x) == 1 % only one event
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
        
        % figure;
        % dendrogram(Z,'ColorThreshold',ops.cutoff)
        % title('Linkage')
        % xlabel('Event')
        % ylabel('Distance')
        % % save figure
        % saveas(gcf, strcat(ops.savedir,sprintf('cluster_ROI%03d',n),'.fig'))
        % saveas(gcf, strcat(ops.savedir,sprintf('cluster_ROI%03d',n),ops.fig_format))
        % close(gcf)

        for i = 1:max(T)
            event_cluster(end+1).ROI = n;
            event_cluster(end).x = peak_x(T==i);
            event_cluster(end).y = peak_y(T==i);
            event_cluster(end).dff_t = peak_t(T==i);
            event_cluster(end).dff = peak_amp(T==i);
            event_cluster(end).n_spikes = nnz(T==i);
            event_cluster(end).ST = 0;

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

    close(gcf)
    assignin("base","event_cluster",event_cluster)
    assignin("base","ROI",ROI)

    %% plot signal for ROI
    figure;
    set(gcf,'Position',[300 300 1500 400])
    title(sprintf('ROI: %03d, Area: %d',n, ROI(n).Area))
    hold on
    plot(ops.t, ROI(n).df)
    scatter(ops.t(ROI(n).t), ROI(n).df(ROI(n).t))
    ylabel('df')
    xlabel('Time [s]')
    % save figure
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),'.fig'))
    saveas(gcf, strcat(ops.savedir_ROIsignal,filesep, sprintf('signal_ROI%03d',n),ops.fig_format))
    close(gcf)
   
end