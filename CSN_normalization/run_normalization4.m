%% init params
clusters_number = 50;
condition_field_name_in_series_matrix = '!Sample_characteristics_name';
first_condition = {'CD4'};
second_condition = {'B Cell'};
first_data_threshold = 4;
second_data_threshold = 4;
FDR_threshold = 0.05;
logFC_threshold = 1;
human_min_iou = 0;
mouse_min_iou = 0;

series_matrix1_filename = 'newGSE107011_series_matrix';
series_matrix2_filename = 'newGSE122597_series_matrix';
exptable1_filename = '4human_pre_norm';
exptable2_filename = '4mouse_pre_norm';
output1_filename = '4human_cmaes';
output2_filename = '4mouse_cmaes';
fig_filename = 'H2-M1_cmaes';


CSN_normalization(clusters_number, condition_field_name_in_series_matrix, first_condition,...
        second_condition, first_data_threshold, second_data_threshold, FDR_threshold, logFC_threshold,...
        human_min_iou, mouse_min_iou, exptable1_filename, exptable2_filename, series_matrix1_filename,...
        series_matrix2_filename, output1_filename, output2_filename, fig_filename);
    
    
human_min_iou = 0.7;
mouse_min_iou = 0.7;
output1_filename = '4human_cmaes_threshold';
output2_filename = '4mouse_cmaes_threshold';
fig_filename = 'H2-M1_cmaes_threshold';


CSN_normalization(clusters_number, condition_field_name_in_series_matrix, first_condition,...
        second_condition, first_data_threshold, second_data_threshold, FDR_threshold, logFC_threshold,...
        human_min_iou, mouse_min_iou, exptable1_filename, exptable2_filename, series_matrix1_filename,...
        series_matrix2_filename, output1_filename, output2_filename, fig_filename);