function ops_multi = config_evoked_multi_spec(ops)
    % total number of images in a multi-spec experiment
    N = 4; % CHANGE THIS TO MATCH YOUR DATA (e.g. 3 for 3-spec data)
    
    % create a cell array to hold ops for each image
    ops_multi = cell(1, N);
    for i = 1:N
        % Initialize each ops structure with default parameters
        ops_multi{i} = default_evoked_spec(ops);
    end

    % Customize parameters for each image as needed
    i = 4; % for the 4th image, 5AP at 20Hz starting at 2s
    ops_multi{i}.first_stim      =  2;                   % Time of first stimulus [s]
    ops_multi{i}.n_stim          =  5;                   % Total number of stimuli
    ops_multi{i}.stim_freq       =  20;                  % Stimulation frequency [Hz]
    ops_multi{i}.len_spike       =  2;                   % Expected response duration [s] from stim onset 


    % After customizing parameters for each image, convert time-based parameters to frame indices
    for i = 1:N
        ops_multi{i} = unit_conversion(ops_multi{i});
    end
    
end

function ops = default_evoked_spec(ops)
    % ========== IMAGE PREPROCESSING OPTIONS ==========
    ops.pre_processing = true;                  % Remove systematic grid line noise (microscope artifact)
    ops.bkg_subtraction = false;                % Background subtraction in time domain
    ops.tophat_max_dff  = true;                 % Tophat filter for uneven illumination correction
    ops.tophat_max_dff_r = 5;                   % Tophat filter radius [pixels]
    ops.use_binary_mask = true;                 % Use user-provided binary mask (ImageJ)
    ops.remove_px_with_no_spikes = true;        % Remove inactive pixels (set false for low SNR data)
    ops.redo_detection = true;                 % Reprocess if results already exist

    % ========== VISUALIZATION & OUTPUT OPTIONS ==========
    ops.plot_pxMap = false;                     % Plot pixel map for each ROI (memory intensive)
    ops.plot_ROI_overall_signal = false;        % Plot averaged signal per ROI (memory intensive)
    ops.plot_ROI_px_signal = false;             % Plot individual pixel signals per ROI (slow)
    ops.fig_format = '.png';                    % Figure format (.png, .pdf, .fig)
    ops.close_fig = true;                       % Close figures after saving (memory efficient)
    ops.visualize = true;                       % Toggle all visualizations (set false for 2-3x speedup)
    ops.save_filtered_tif = false;              % Save ImJFig1_Max_dFoF_TopHatFiltered.tif (tophat filter mask)
    ops.save_px_mask = false;                   % Save ImJFig2_PxMask.tif (pixel mask)

    % ========== TEMPORAL FILTERING PARAMETERS ==========
    ops.baseline_percentage = 0.2;              % Baseline percentile for drift estimation [0-1]
    ops.sl_window           = 3;                % Sliding window for baseline [s]
    ops.sl_window_ST        = 5;                % Sliding window for spike train [s]

    % ========== PIXEL-LEVEL FILTERING THRESHOLDS ==========
    ops.filter_by_slope = true;                 % Enable SNR-based pixel filtering
    ops.m_mode    = 'SD';                       % Mode: 'SD' (standard deviation) or 'absolute'
    ops.m_thres   = 3;                          % Threshold multiplier (3 = 3x std above baseline)
    ops.SNR_thres = 5;                          % Signal-to-noise ratio threshold for pixel inclusion

    % ========== BASELINE DRIFT DETECTION ==========
    ops.BL_drift_thres = 0.05;                  % Max RMS error for linear baseline fit (remove if exceeded)

    % ========== SPONTANEOUS EVENT DETECTION PARAMETERS ==========
    % These thresholds define event characteristics for spontaneous activity
    ops.spontaneous.rising_time_thres  = 0.05; % Minimum event rise time [s] - filters slow drifts
    ops.spontaneous.maxISI             = 0.08; % Maximum inter-spike interval [s] - fusion threshold
    ops.spontaneous.findpeak_window    = 0.07; % Time window for peak detection [s]
    ops.spontaneous.dfof_MinPeakHeight = 3;    % Minimum peak height [multiples of pixel-level std]
    ops.spontaneous.peakWidth          = 0.4;  % Expected event width [s] in ΔF/F signal
    ops.spontaneous.MinPeakWidth       = 0.02; % Minimum peak width for 1D detection [s]
    ops.spontaneous.MinPeakHeight      = 150;  % Minimum peak height for 2D spatial detection

    % ========== EXPERIMENT TYPE & STIMULATION PARAMETERS ==========
    ops.experiment_type = "evoked";             % Type: "spontaneous" or "evoked"
    % NOTE: evoked parameters below are unused for spontaneous experiments % 15AP at 1Hz starting at 3s
    ops.first_stim      =  3;                   % Time of first stimulus [s]
    ops.n_stim          = 15;                   % Total number of stimuli
    ops.stim_freq       =  1;                   % Stimulation frequency [Hz]
    ops.len_spike       =  2;                   % Expected response duration [s] from stim onset

    % ========== EVOKED EVENT DETECTION PARAMETERS (if experiment_type == "evoked") ==========
    ops.evoked.rising_time_thres  = 0.15;      % Minimum response rise time [s]
    ops.evoked.maxISI             = 0.08;      % Maximum inter-spike interval [s]
    ops.evoked.findpeak_window    = 0.15;      % Time window for peak search [s]
    ops.evoked.dfof_MinPeakHeight = 2.5;       % Minimum response amplitude [multiples of std]
    ops.evoked.peakWidth          = 2;         % Expected response width [s]
    ops.evoked.MinPeakWidth       = 0.02;      % Minimum peak width for 1D detection [s]
    ops.evoked.MinPeakHeight      = 0.1;       % Minimum peak height for 2D detection

    % ========== ACQUISITION & ANALYSIS PARAMETERS ==========
    ops.fs = 100;                               % Frame rate [Hz] - MUST match data acquisition rate

    % ROI SIZE FILTERING (optional - comment out to disable)
    ops.Area_thres_max = 50;                    % Maximum ROI area [pixels] - filters noise clusters
    % ops.Area_thres_min = 1;                  % Minimum ROI area [pixels] - uncomment if needed

    % ========== SPATIAL CLUSTERING PARAMETERS ==========
    ops.cutoff = 5;                             % Distance cutoff for clustering [pixels]
                                                % ROIs within this distance + synchronized timing = merged
    ops.timing_tolerance = 2;                   % Timing tolerance for spike synchronization [frames]
                                                % ROI spikes must align within this window to be merged

    % ========== SPIKE TRAIN ANALYSIS PARAMETERS ==========
    ops.ST.option           = "findpeaks";      % Method: "MLspike" (recommended) or "findpeaks"
    ops.ST.MinPeakHeight    = 0.05;             % Minimum spike amplitude
    ops.ST.sumOfPeak_window = 0.5;              % Window for spike summation [s]
    ops.ST.gaussian_window  = 0.1;              % Gaussian smoothing window [s]
    ops.ST.gap_thres        = 0.3;              % Gap threshold for spike identification [s]
    ops = spike_train_par(ops);

    
end

function ops = unit_conversion(ops)
    % Unit conversion: Convert time-based parameters from seconds to frames
    % This allows all operations to work with frame indices instead of time values
    ops.sl_window          = round(ops.sl_window * ops.fs);              % Baseline window [frames]
    ops.sl_window_ST       = round(ops.sl_window_ST * ops.fs);           % Spike train window [frames]

    % Convert spontaneous event detection thresholds to frames
    ops.spontaneous.rising_time_thres  = round(ops.spontaneous.rising_time_thres * ops.fs);
    ops.spontaneous.maxISI             = round(ops.spontaneous.maxISI * ops.fs);
    ops.spontaneous.findpeak_window    = round(ops.spontaneous.findpeak_window * ops.fs);
    ops.spontaneous.peakWidth          = round(ops.spontaneous.peakWidth * ops.fs);
    ops.spontaneous.MinPeakWidth       = round(ops.spontaneous.MinPeakWidth * ops.fs);

    % Convert evoked parameters if running evoked experiment
    if ops.experiment_type == "evoked"
        ops.(ops.experiment_type).rising_time_thres  = round(ops.(ops.experiment_type).rising_time_thres * ops.fs);
        ops.(ops.experiment_type).maxISI             = round(ops.(ops.experiment_type).maxISI * ops.fs);
        ops.(ops.experiment_type).findpeak_window    = round(ops.(ops.experiment_type).findpeak_window * ops.fs);
        ops.(ops.experiment_type).peakWidth          = round(ops.(ops.experiment_type).peakWidth * ops.fs);
        ops.(ops.experiment_type).MinPeakWidth       = round(ops.(ops.experiment_type).MinPeakWidth * ops.fs);

        % Calculate stimulus timing in frames
        ops.stim_time = (0:ops.n_stim-1)/ops.stim_freq + ops.first_stim;
        ops.stim_frames = round(ops.stim_time * ops.fs);
        ops.stim_pk_search_range = arrayfun(@(x) x:x+ops.(ops.experiment_type).findpeak_window, ops.stim_frames, 'UniformOutput', false);
        ops.stim_pk_search_range = reshape(ops.stim_pk_search_range,[],1);
        ops.stim_pk_search_range = cell2mat(ops.stim_pk_search_range);
        ops.len_spike = round(ops.len_spike * ops.fs);
    end

    % Convert spike train parameters to frames
    ops.ST.sumOfPeak_window = round(ops.ST.sumOfPeak_window * ops.fs);
    ops.ST.gaussian_window  = round(ops.ST.gaussian_window * ops.fs);
    ops.ST.gap_thres        = round(ops.ST.gap_thres * ops.fs);
end












