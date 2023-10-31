%% Sliding window filter
function baseline_pc_median = sliding_window_filter(data, baseline_percentage, window)
    tic;
    disp('Applying sliding window filter...')
    zeropad_len = window;
    data = zero_padding(data, zeropad_len);
    T = size(data,1);
    baseline_pc_median = zeros(size(data));
    for k=1:T %for all timepoints
        kymo_sample = data(max(1,k-window/2):min(T,k+window/2),:); %take window around each timepoint
        sortedVals = sort(kymo_sample); %sort in ascending order
        baseline_sorted_percentage=sortedVals(1:round(baseline_percentage*size(kymo_sample,1)),:); %take a specified percentage of these values
        baseline_pc_median(k,:) = median(baseline_sorted_percentage); %find median of this percentage
    end
    baseline_pc_median = remv_zero_padding(baseline_pc_median, zeropad_len);
    toc
end
