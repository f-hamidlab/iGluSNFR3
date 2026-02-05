function ops = config_spontaneous(ops)
    % ========== IMAGE PREPROCESSING OPTIONS ==========
    ops.pre_processing = true;                  % Remove systematic grid line noise (microscope artifact)
    ops.bkg_subtraction = false;                % Background subtraction in time domain
    ops.tophat_max_dff  = true;                 % Tophat filter for uneven illumination correction
    ops.tophat_max_dff_r = 12;                  % Tophat filter radius [pixels]
    ops.use_binary_mask = false;                % Use user-provided binary mask (ImageJ)
    ops.remove_px_with_no_spikes = true;        % Remove inactive pixels (set false for low SNR data)
    ops.redo_detection = true;                  % Reprocess if results already exist

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
    ops.m_thres   = 4;                          % Threshold multiplier (4 = 4x std above baseline)
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
    ops.spontaneous.MinPeakHeight      = 0.01; % Minimum peak height for 2D spatial detection

    % ========== EXPERIMENT TYPE & STIMULATION PARAMETERS ==========
    ops.experiment_type = "spontaneous";       % Type: "spontaneous" or "evoked"
    % NOTE: evoked parameters below are unused for spontaneous experiments
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
    ops.timing_tolerance = 1;                   % Timing tolerance for spike synchronization [frames]
                                                % ROI spikes must align within this window to be merged

    % ========== SPIKE TRAIN ANALYSIS PARAMETERS ==========
    ops.ST.option           = "MLspike";       % Method: "MLspike" (recommended) or "findpeaks"
    ops.ST.MinPeakHeight    = 0.05;             % Minimum spike amplitude
    ops.ST.sumOfPeak_window = 0.5;              % Window for spike summation [s]
    ops.ST.gaussian_window  = 0.1;              % Gaussian smoothing window [s]
    ops.ST.gap_thres        = 0.3;              % Gap threshold for spike identification [s]
    ops = spike_train_par(ops);

    
end