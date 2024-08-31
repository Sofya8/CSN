function [human_weight_for_human, mouse_weight_for_human,...
    mouse_weight_for_mouse, human_weight_for_mouse,...
    human_alphas, human_betas, mouse_alphas, mouse_betas] = vec_to_params(x_vec, clusters_number)

	human_weight_for_human = x_vec(1)+1;
    mouse_weight_for_human = x_vec(2);
    mouse_weight_for_mouse = x_vec(3)+1;
    human_weight_for_mouse = x_vec(4);
    human_alphas = x_vec(5:clusters_number+4)+1; 
    human_betas = x_vec(clusters_number+5:clusters_number*2+4);
    mouse_alphas = x_vec(clusters_number*2+5:clusters_number*3+4)+1; 
    mouse_betas = x_vec(clusters_number*3+5:end);
end