function CSN_normalization(clusters_number, condition_field_name_in_series_matrix, first_condition,...
    second_condition, first_data_threshold, second_data_threshold, FDR_threshold, logFC_threshold,...
    human_min_iou, mouse_min_iou, exptable1_filename, exptable2_filename, series_matrix1_filename,...
    series_matrix2_filename, output1_filename, output2_filename, fig_filename)

    parameters_num = clusters_number*4+4;
    condition_names = [first_condition, second_condition];
    fig_filename = strcat('output/',fig_filename);

    series_matrix1_filename = strcat('data/',series_matrix1_filename);
    series_matrix2_filename = strcat('data/',series_matrix2_filename);
    fstSeriesTable = readcell(series_matrix1_filename);
    secSeriesTable = readcell(series_matrix2_filename);
    % 

    exptable1_filename = strcat('data/',exptable1_filename);
    exptable2_filename = strcat('data/',exptable2_filename);
    fstExptable = readtable(exptable1_filename);
    secExptable = readtable(exptable2_filename);
    origin_data_human = table2array(fstExptable(:,2:end));
    origin_data_mouse = table2array(secExptable(:,2:end));



    % % 
    [samples_name_human_cond1, samples_name_human_cond2,...
        samples_name_mouse_cond1, samples_name_mouse_cond2] = ...
                        extract_samples_and_build_group_tables(fstSeriesTable, secSeriesTable,...
                                                               condition_field_name_in_series_matrix,...
                                                               condition_names);

    HC1_table = [fstExptable(:,1) fstExptable(:,samples_name_human_cond1)];
    HC2_table = [fstExptable(:,1) fstExptable(:,samples_name_human_cond2)];
    MC1_table = [secExptable(:,1) secExptable(:,samples_name_mouse_cond1)];
    MC2_table = [secExptable(:,1) secExptable(:,samples_name_mouse_cond2)];
    HC1_size = size(HC1_table,2)-1;
    HC2_size = size(HC2_table,2)-1;
    MC1_size = size(MC1_table,2)-1;
    MC2_size = size(MC2_table,2)-1;
 
    origin_DEGs_HC1_HC2 = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               HC2_table, HC2_size, first_data_threshold, FDR_threshold, logFC_threshold);
	origin_DEGs_MC1_MC2 = find_DEGs(MC1_table, MC1_size, second_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold); 
    
	merged_data = [origin_data_human, origin_data_mouse];
    index = run_fuzzy_k_means(merged_data, clusters_number);


    normalization_function = @(x_vec)loss_func(x_vec, fstExptable, secExptable, origin_data_human,...
                                               origin_data_mouse,...
                                               index, clusters_number, samples_name_human_cond1,...
                                               samples_name_human_cond2, samples_name_mouse_cond1,...
                                               samples_name_mouse_cond2,  first_data_threshold,...
                                               second_data_threshold, FDR_threshold, logFC_threshold,...
                                               origin_DEGs_HC1_HC2, origin_DEGs_MC1_MC2, human_min_iou, mouse_min_iou);

   cmaes_res = cmaes(normalization_function, parameters_num, -2, 2, 50, fig_filename);
   
    normalization_function(cmaes_res);
	[new_data_human, new_data_mouse] = transform_datasets_accord_to_vec(origin_data_human,...
                                                                        origin_data_mouse,...
                                                                        index, cmaes_res, clusters_number);
    
    fstExptable(:,2:end) = array2table(new_data_human);
    secExptable(:,2:end) = array2table(new_data_mouse);

    writetable(fstExptable, strcat('output/',output1_filename));
    writetable(secExptable, strcat('output/',output2_filename));
end
function cellArray = readcell(path)
    T = readtable(path);
    cellArray = table2cell(T);
end