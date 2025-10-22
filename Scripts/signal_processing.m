function [signal_baseline, signal_df, signal_dfof, signal_dfof_movemean, signal_edge] = signal_processing(signal_raw, ops, option)
    switch option
        case "spontaneous"
            signal_baseline = sliding_window_filter(signal_raw, ops.baseline_percentage, ops.sl_window);
            signal_df = signal_raw - signal_baseline;
            signal_baseline = signal_baseline + median(signal_df,1);
            signal_df = signal_df - median(signal_df,1);
            signal_dfof = signal_df./signal_baseline;
            signal_edge = Gaussian_Derivative_Filter_padded(signal_dfof, ops.m_mode, ops.m_thres);
            signal_dfof_movemean = movmean(signal_dfof, 5, 1); % window of 5
            
        case "evoked"
            signal_baseline = sliding_window_filter(signal_raw, ops.baseline_percentage, ops.sl_window);
            signal_df = signal_raw - signal_baseline;
            signal_baseline = signal_baseline + median(signal_df,1);
            signal_df = signal_df - median(signal_df,1);
            signal_dfof = signal_df./signal_baseline;
            signal_edge = Gaussian_Derivative_Filter_padded(signal_dfof, ops.m_mode, ops.m_thres);
            signal_dfof_movemean = movmean(signal_dfof, 5, 1); % window of 5

        case "evoked_mask"
            f1 = ops.stim_frames(1)-1;
            f2 = ops.stim_frames(end)+ops.len_spike+1;
            signal_baseline = zeros(size(signal_raw));
            signal_baseline(1:f1) = prctile(signal_raw(1:f1),ops.baseline_percentage*100);
            signal_baseline(f2:end) = prctile(signal_raw(f2:end),ops.baseline_percentage*100);
            % interpolation for frames inbetween
            signal_baseline(f1:f2) = interp1([f1 f2], [signal_baseline(f1) signal_baseline(f2)], f1:f2);
            
            signal_df = signal_raw - signal_baseline;
            signal_dfof = signal_df./signal_baseline;
            signal_edge = Gaussian_Derivative_Filter_padded(signal_dfof, ops.m_mode, 1); % assume noise is small for overall trend. Hence threshold set to 1.
            signal_dfof_movemean = movmean(signal_dfof, 5, 1); % window of 5

    end
    
end
%%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%                                                                         %
%                            local functions                              %
%                                                                         %       
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%% Gaussian filter
function edge = Gaussian_Derivative_Filter_padded(data, mode, thres)
    tic;
    disp('Applying 1st order Gaussian filter...')

    % sigma for guassian filter
    sigma_y = 2;
    
    % length of zero-padding
    zeropad_len = 50;
    
    % zero-padding
    edge = zero_padding(data, zeropad_len);
    
    % 1st Order Derivative 1D Gaussian filter, detects slopes in time dimension
    edge = Gaussian_Derivative_Filter(edge, sigma_y);
    
    % remove zero-padding
    edge = remv_zero_padding(edge, zeropad_len);
    
    % correct for direction of the slope
    edge = -edge;
    
    % thresholding to remove slow flucturations in signal
    switch mode
        case 'absolute'
            edge(edge<thres)=0;
        case 'SD'
            STD = std(edge);
            edge(edge<thres*STD)=0;
    end

    toc
end

function im = Gaussian_Derivative_Filter(im, sigma_y)
    y_values = -ceil(4*sigma_y):ceil(4*sigma_y);
    filter_y = calc1stOrderDerivative1DGaussian(y_values, sigma_y)';
    im = imfilter(im, filter_y,'symmetric');
end

function [ D_values ] = calc1stOrderDerivative1DGaussian(x_values, sigma)
    % Getting the 1D Gaussian values
    G_values = calc1DGaussian(x_values, sigma);
    
    % Obtaining the 1st order derivative of the 
    D_values = -2.*x_values./2./sigma^2.*G_values;
end

function [ G_values ] = calc1DGaussian(x_values, sigma)
    G_values = 1/sigma/sqrt(2*pi).*exp(-x_values.^2./2./sigma^2);
end
