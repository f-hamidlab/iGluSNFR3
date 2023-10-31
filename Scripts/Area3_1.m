% adjust spike train detection parameters (optional)
ops = spike_train_par(ops);

% add spike trains
spike_train(67); % not sure
spike_train(70); % not sure
spike_train(124);
spike_train(125);

% remove ROIs with only noise
event_cluster = remove_ROI(event_cluster, 9);
event_cluster = remove_ROI(event_cluster, 10);
event_cluster = remove_ROI(event_cluster, 40);
event_cluster = remove_ROI(event_cluster, 41);
% event_cluster = remove_ROI(event_cluster, 71); % not sure
event_cluster = remove_ROI(event_cluster, 87);

% add events
add_event(31, 1.94)

% remove events (noise in an ROI that you want to keep)
% remove_event(31, 1.94) % (example only)

%% event statistics
stats = event_stats(event_cluster, ops);

% show label mask
show_label_mask(event_cluster, ROI, ops)

% save data as .mat
save_data()
