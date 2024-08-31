function values_to_minimize = loss_func(x_vec, human_table, mouse_table, data_human, data_mouse,...
                                                       index, clusters_number, samples_name_human_cond1,...
                                                       samples_name_human_cond2, samples_name_mouse_cond1,...
                                                       samples_name_mouse_cond2,  first_data_threshold,...
                                                       second_data_threshold, FDR_threshold, logFC_threshold,...
                                                       human_origin_DEGS, mouse_origin_DEGS, human_min_iou, mouse_min_iou)

    [new_data_human, new_data_mouse] = transform_datasets_accord_to_vec(data_human, data_mouse,...
                                                                        index, x_vec, clusters_number);
     human_table(:, 2:end) = array2table(new_data_human);
     mouse_table(:, 2:end) = array2table(new_data_mouse);
    
    HC1_table = [human_table(:,1), human_table(:,samples_name_human_cond1)];
    HC2_table = [human_table(:,1), human_table(:,samples_name_human_cond2)];
    MC1_table = [mouse_table(:,1), mouse_table(:,samples_name_mouse_cond1)];
    MC2_table = [mouse_table(:,1), mouse_table(:,samples_name_mouse_cond2)];
    HC1_size = size(HC1_table,2)-1;
    HC2_size = size(HC2_table,2)-1;
    MC1_size = size(MC1_table,2)-1;
    MC2_size = size(MC2_table,2)-1;
 
    HC1_HC2 = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               HC2_table, HC2_size, first_data_threshold, FDR_threshold, logFC_threshold);
	MC1_MC2 = find_DEGs(MC1_table, MC1_size, second_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold);
	HC1_MC2 = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold);
	MC1_HC2 = find_DEGs(MC1_table, MC1_size, second_data_threshold,....
                            HC2_table, HC2_size, first_data_threshold, FDR_threshold, logFC_threshold);
    HC1_MC1 = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               MC1_table, MC1_size, second_data_threshold, FDR_threshold, logFC_threshold);
	HC2_MC2 = find_DEGs(HC2_table, HC2_size, first_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold);
                           
	list_of_tables = {HC1_HC2, MC1_MC2, HC1_MC2, MC1_HC2, HC1_MC1, HC2_MC2};
    
    number_of_lists = length(list_of_tables);
    precentMatrix = zeros(number_of_lists,number_of_lists);
    for i = 1 : number_of_lists
       precentMatrix(i, i) = 1;
       for j = i+1 : number_of_lists
          [precent_genes_fst_table_in_intrsct,...
                  precent_genes_sec_table_in_intrsct] = intersect_between_gene_lists(list_of_tables{i},...
                                                            list_of_tables{j});
          precentMatrix(i, j) = precent_genes_fst_table_in_intrsct;
          precentMatrix(j, i) = precent_genes_sec_table_in_intrsct;
       end
    end
    precentMatrix(isnan(precentMatrix))= 0;
    
    square_value = (sum(sum(precentMatrix(1:4,1:4)))-4)/12;
    frame_value = (sum(sum(precentMatrix(1:4,5:6))) + sum(sum(precentMatrix(5:6,1:4))))/16;
    CSC_index = square_value/frame_value;

    
    IOU_human = intersection_over_union(list_of_tables{1}, human_origin_DEGS);
    IOU_mouse = intersection_over_union(list_of_tables{2}, mouse_origin_DEGS);
    
     values_to_minimize = 1000 * (max((human_min_iou - IOU_human),0) + max((mouse_min_iou - IOU_mouse),0)) - CSC_index - IOU_human - IOU_mouse;

end