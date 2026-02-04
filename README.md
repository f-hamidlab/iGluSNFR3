# iGluSNFR3
Detection of Synaptic Activities in Widefield Microscopy

For installation and details of the toolbox, please refer to the [wiki page](https://github.com/f-hamidlab/iGluSNFR3/wiki/02-Using-the-Scripts).

## Designed for:
- single plane image 
- [iGluSNFR3 probe] (https://www.nature.com/articles/s41592-023-01863-6)
- multi-cellular
- spontaneous synaptic activities

## Dependencies:
1) MLspike toolbox (https://github.com/MLspike)
2) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html, https://bio-formats.readthedocs.io/en/v7.0.0/developers/matlab-dev.html)

## Main Pipeline Scripts:
- multi_cell_activity_detection_pipeline.m - Main pipeline for synaptic activity detection
- multi_cell_activity_detection_pipeline_spontaneous.m - Spontaneous activity detection
- multi_cell_activity_detection_pipeline_evoked.m - Evoked activity detection
- multi_cell_activity_detection_pipeline_evoked_matching.m - Evoked activity detection with cluster matching

## Signal Processing & Analysis Scripts:
- signal_processing.m - Core signal processing routines
- findpeaks_2d.m - 2D peak detection
- sliding_window_filter.m - Sliding window filtering
- cluster_events.m - Event clustering
- knee_pt.m - Knee point detection
- px_properties.m - Pixel properties analysis
- mask_properties_evoked_overall.m - Mask properties for evoked responses

## ROI & Spatial Analysis:
- ROI_pxMap_allTime.m - ROI pixel mapping across time
- Area3_1.m - Area analysis
- matching_clusters.m - Cluster matching across conditions
- matching_clusters_with_image4.m - Cluster matching with image registration

## Spike Train & Event Analysis:
- spike_train.m - Spike train generation and analysis
- spike_train_par.m - Parallel spike train analysis
- ST_function.m - Spike train functions
- event_stats.m - Event statistics calculation
- plotSpikeRaster.m - Spike raster plotting
- rasterplot.m - Raster plot visualization

## Data I/O & Visualization:
- save_data.m - Data saving utilities
- saveastiff.m - TIFF file saving
- show_label_mask.m - Labeled mask visualization
- show_label_mask_with_text.m - Labeled mask with text labels
- plot_signal_pixel.m - Pixel signal plotting
- loadBioFormats.m - Bio-Formats image loading

## Data Curation & Editing:
- add_event.m - Add events manually
- remove_event.m - Remove events
- remove_ROI.m - Remove ROI regions

