function mask = mask_properties_evoked_overall(mask, ops)
    [pks,locs] = findpeaks(mask.signal_edge);
    % skip the first and last 5 frames
    pks = pks(locs>5);
    locs = locs(locs>5);
    pks = pks(locs<ops.Nt-5);
    locs = locs(locs<ops.Nt-5);

    mask.slope_pk = pks;
    mask.slope_t = locs;

    % ignore spikes during calculation of signal STD
    px_dfof_std = mask.signal_dfof;
    
    for j=1:length(locs)
        y=locs(j);
        px_dfof_std(y:y+ops.(ops.experiment_type).peakWidth) = NaN;

    end
    mask.STD = std(px_dfof_std(:),'omitnan');

    for j = 1:length(locs) % for each frame where peak slope is located
        t_start = locs(j);
        t_end = min([ops.Nt, t_start+ops.(ops.experiment_type).findpeak_window]);

        [~, dff_t] = findpeaks(mask.signal_dfof(t_start:t_end),'NPeaks',1, ...
            'MinPeakHeight', ops.(ops.experiment_type).dfof_MinPeakHeight*mask.STD, ...
            'WidthReference','halfheight',...
            'MinPeakWidth',ops.(ops.experiment_type).MinPeakWidth);
        
        % find the first zero-crossing
        % if the first zero-crossing is over threshold, remove pt
        if ~isempty(dff_t)
            dff_t = dff_t+t_start-1;
            
            if ~isempty(mask.dff_t)
                last_pk = mask.dff_t(end);
                dis = dff_t-last_pk;
                if dis<ops.(ops.experiment_type).maxISI
                    mask.df(end+1) = mask.signal_df(dff_t);
                    mask.dff_t(end+1) = dff_t;
                    mask.dff(end+1) = mask.signal_dfof(dff_t);
                    continue
                else 
                    % check time of rising slope
                    ipt = find(mask.signal_edge(1:t_start)<=0,1,"last");
                    if ~isempty(ipt)
                        mask.ipt(end+1) = ipt;
                    end
                    % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                    rising_time = dff_t-ipt;
                end
            else 
                % check time of rising slope
                ipt = find(mask.signal_edge(1:t_start)<=0,1,"last");
                if ~isempty(ipt)
                    mask.ipt(end+1) = ipt;
                end
                % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                rising_time = dff_t-ipt;
            end
            
            if rising_time>ops.(ops.experiment_type).rising_time_thres
                continue
            end
            
            mask.dff_t(end+1) = dff_t;
            mask.df(end+1) = mask.signal_df(dff_t);
            mask.dff(end+1) = mask.signal_dfof(dff_t);

            mask.dff_from_ipt(end+1) = (mask.signal_raw(dff_t) - mask.signal_raw(ipt)) / mask.signal_baseline(dff_t);

        end
    end
    mask.n_spikes = length(mask.dff_t);
end