function new_fst_data = standrdize(origin_fst_data, origin_sec_data,...
                                     fst_dataset_weight, sec_dataset_weight)
    % pre - processing
    human_mean = mean(origin_fst_data, 2);
    human_std = std(origin_fst_data')';
    mouse_mean = mean(origin_sec_data, 2);
    mouse_std = std(origin_sec_data')';
    
    human_mouse_mean = (human_mean*fst_dataset_weight + mouse_mean*sec_dataset_weight)/(fst_dataset_weight + sec_dataset_weight);
    human_mouse_std = (human_std*fst_dataset_weight + mouse_std*sec_dataset_weight)/(fst_dataset_weight + sec_dataset_weight);

    new_fst_data = (origin_fst_data - human_mean)./human_std.*human_mouse_std+human_mouse_mean;
    new_fst_data(isnan(new_fst_data)) = 0;
end