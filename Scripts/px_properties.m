%% px_properties
% Calculates event detection properties for individual pixels
%
% DESCRIPTION:
%   Analyzes temporal dynamics of individual pixels including peak slopes,
%   signal variability, and event timing. Identifies peak times and amplitudes
%   using derivative-based peak detection with noise estimation.
%
% USAGE:
%   1) px = px_properties(px, ops, px_edge, px_df, px_dfof, px_dfof_movemean)
%
% INPUTS:
%   - px: (struct) pixel data structure
%   - ops: (struct) options and parameters
%   - px_edge: (numeric) temporal derivative of pixel signal (edge detection)
%   - px_df: (numeric) delta F signal for pixel
%   - px_dfof: (numeric) delta F/F signal for pixel
%   - px_dfof_movemean: (numeric) moving average of delta F/F
%
% OUTPUTS:
%   - px: (struct) updated pixel structure with properties:
%       .slope_pk: peak slope magnitudes
%       .slope_t: times of peak slopes
%       .STD: signal noise level (std of non-spike regions)
%       .dff_t: times of detected events
%       .dff: amplitudes of detected events
%       .df: delta F amplitudes
%       .SNR: signal-to-noise ratio
%
% Last updated: 2026-02-03 15:30

function px = px_properties(px, ops, px_edge, px_df, px_dfof, px_dfof_movemean)
    [pks,locs] = findpeaks(px_edge);
    % skip the first and last 5 frames
    pks = pks(locs>5);
    locs = locs(locs>5);
    pks = pks(locs<ops.Nt-5);
    locs = locs(locs<ops.Nt-5);

    px.slope_pk = pks;
    px.slope_t = locs;

    % ignore spikes during calculation of signal STD
    px_dfof_std = px_dfof;
    
    for j=1:length(locs)
        y=locs(j);
        px_dfof_std(y:y+ops.(ops.experiment_type).peakWidth) = NaN;

    end
    px.STD = std(px_dfof_std(:),'omitnan');

    for j = 1:length(locs) % for each frame where peak slope is located
        t_start = locs(j);
        t_end = min([ops.Nt, t_start+ops.(ops.experiment_type).findpeak_window]);

        [~, dff_t] = findpeaks(px_dfof(t_start:t_end),'NPeaks',1, ...
            'MinPeakHeight', ops.(ops.experiment_type).dfof_MinPeakHeight*px.STD, ...
            'WidthReference','halfheight',...
            'MinPeakWidth',ops.(ops.experiment_type).MinPeakWidth);
        
        % find the first zero-crossing
        % if the first zero-crossing is over threshold, remove pt
        if ~isempty(dff_t)
            dff_t = dff_t+t_start-1;
            
            if ~isempty(px.dff_t)
                last_pk = px.dff_t(end);
                dis = dff_t-last_pk;
                if dis<ops.(ops.experiment_type).maxISI
                    px.df(end+1) = px_df(dff_t);
                    px.dff_t(end+1) = dff_t;
                    px.dff(end+1) = px_dfof(dff_t);
                    continue
                else 
                    % check time of rising slope
                    STD_crossing = find(px_dfof_movemean(1:dff_t)<3*px.STD,1,"last");
                    t_find_ch_pt = max([STD_crossing-50, round((last_pk+dff_t)/2), 1]);
                    ipt = findchangepts(px_dfof_movemean(t_find_ch_pt:dff_t),Statistic="linear", MaxNumChanges=1);
                    ipt = ipt + t_find_ch_pt;
                    if ~isempty(ipt)
                        px.ipt(end+1) = ipt;
                    end
                    % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                    rising_time = dff_t-ipt;
                end
            else 
                % check time of rising slope
                STD_crossing = find(px_dfof_movemean(1:dff_t)<3*px.STD,1,"last");
                ipt = findchangepts(px_dfof_movemean(max([STD_crossing-50, 1]):dff_t),Statistic="linear", MaxNumChanges=1);
                ipt = ipt + max([STD_crossing-50, 1])-1;
                if ~isempty(ipt)
                    px.ipt(end+1) = ipt;
                end
                % zero_crossing = find(signal_dfof(1:dff_locs,i)<0,1,"last");
                rising_time = dff_t-ipt;
            end
            
            if rising_time>ops.(ops.experiment_type).rising_time_thres
                continue
            end

            px.df(end+1) = px_df(dff_t);
            px.dff_t(end+1) = dff_t;
            px.dff(end+1) = px_dfof(dff_t);

        end
    end
    px.n_spikes = length(px.dff_t);
end