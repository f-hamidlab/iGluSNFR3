%% ST_function
% Applies spike train classification threshold to ROI events
%
% DESCRIPTION:
%   Determines which ROIs contain true spike trains versus noise based on
%   spike width and peak count characteristics. Applies a linear threshold
%   to classify ROIs as having valid spike trains.
%
% USAGE:
%   1) ROI = ST_function(ROI, ops)
%
% INPUTS:
%   - ROI: (struct) ROI data as structure array
%   - ops: (struct) options and parameters including .ST field with spike train settings
%
% OUTPUTS:
%   - ROI: (struct) updated ROI structure with added fields:
%       .ST_width: spike train width measure
%       .ST_NoP: number of peaks in spike train window
%       .ST: (logical) classification as true spike train
%
% NOTES:
%   Prints total count of ROIs classified as having spike trains
%
% Last updated: 2026-02-03 15:30

function ROI = ST_function(ROI, ops)
    
    P = arrayfun(@(x) getNoiseLevel(x.df, x.t, ops.ST),ROI,'UniformOutput',false);  
    [ROI(:).ST_width] = P{:};

    all_peaks_num = zeros(size(ROI,1),1); % Determine number of peaks in spike train
    for i = 1:length(ROI)
        sig = zeros (1,ops.Nt);
        sig(ROI(i).t) = 1;
        a = movsum(sig, ops.ST.sumOfPeak_window);
        peak_num = (max(a)); % Number of peaks within the window 
        all_peaks_num(i) = peak_num;
    end 
    M = all_peaks_num;
    M = num2cell(M);
    [ROI(:).ST_NoP] = M{:};

    x2 = [ROI(:).ST_width];
    y2 = [ROI(:).ST_NoP];
    ST = (-0.045*x2)+5 < y2; % Spike train threshold 
    ST = num2cell(ST);
    [ROI(:).ST] = ST{:};
    total_st = sum([ROI(:).ST]); % Total number of spike trains

    formatSpec = 'Number of ROIs with spike trians = %d \n';
    fprintf(formatSpec, total_st); % Display total number of spike trains 

end 

function width = getNoiseLevel(noise_signal, t, ST) % Determine width of spike train 
    % Defining peaks 
    smoothed_signal = smoothdata(noise_signal,"gaussian", ST.gaussian_window);
    amplitude = smoothed_signal(t); 
    a = amplitude / 2; % Midpoint of peak
    all_signal = zeros(1,length(smoothed_signal)); 
    for i = 1 : length(t)
        x_start = find(smoothed_signal(1:t(i)) < a(i),1,'last'); % Start point of peak
        x_end = find(smoothed_signal(t(i):end) < a(i),1,'first');
        x_end = x_end+t(i)-1; % End point of peak 
        all_signal(x_start:x_end) = 1; % All points within peak = 1
    end
    % Determiningthe gap between peaks
    f = find(diff([0,~all_signal,0] == 1)); 
    s = f(1:2:end-1); % Start of the gap
    l = f(2:2:end)-s; % End of the gap 

    % If the gap between 2 peaks is < 30 = spike train 
    for i = 1:length(l)
        if  l(i) < ST.gap_thres
            all_signal(s(i):(s(i)+l(i))) = 1;
        end 
    end 
    f = find(diff([0,all_signal,0] == 1));
    s = f(1:2:end-1);  
    y = f(2:2:end)-s; 
    width = max(y); % Establishing widest spike train 
end
