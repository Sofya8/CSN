function [new_data_human, new_data_mouse] = transform_datasets_accord_to_vec(data_human, data_mouse,...
                                                                             index, x_vec, clusters_number)

    [human_weight_for_human, mouse_weight_for_human,...
        mouse_weight_for_mouse, human_weight_for_mouse,...
        human_alphas, human_betas, mouse_alphas, mouse_betas] = vec_to_params(x_vec, clusters_number);
    
    standardized_data_human = standrdize(data_human, data_mouse, human_weight_for_human, mouse_weight_for_human);
    new_data_human = move_blocks(standardized_data_human, index, human_alphas, human_betas, clusters_number);
	
    standardized_data_mouse = standrdize(data_mouse, data_human, mouse_weight_for_mouse, human_weight_for_mouse);
    new_data_mouse = move_blocks(standardized_data_mouse, index, mouse_alphas, mouse_betas, clusters_number);
end