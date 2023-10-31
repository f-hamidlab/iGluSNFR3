# iGluSNFR3
Detection of Synaptic Activities in 2-Photon Microscopy

For installation and details of the toolbox, please refer to the [wiki page](https://github.com/f-hamidlab/iGluSNFR3/wiki/02-Using-the-Scripts).

## Designed for:
- single plane image 
- iGluSNFR3 probe
- multi-cellular
- spontaneous synaptic activities

## Dependencies:
1) MLspike toolbox (https://github.com/MLspike)
2) Bio-Format Toolbox (https://bio-formats.readthedocs.io/en/v7.0.0/users/matlab/index.html, https://bio-formats.readthedocs.io/en/v7.0.0/developers/matlab-dev.html)
3) Additional scripts:
      - event_stats.m
      - findpeaks_2d.m
      - remove_ROI.m
      - remv_zero_padding.m
      - save_data.m
      - show_label_mask.m
      - sliding_window_filter.m
      - zero_padding.m

## Additional scripts/functions useful for manual curation of data:
- add_event.m
- event_stats.m
- plot_signal_pixel.m
- remove_event.m
- remove_ROI.m
- ROI_pxMap_allTime.m
- save_data.m
- show_label_mask.m
- spike_train.m
- spike_train_par.m

