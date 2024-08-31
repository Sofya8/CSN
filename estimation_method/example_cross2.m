% clc
% close all

condition_field_name_in_series_matrix = '!Sample_characteristics_name';
first_condition = {'CD4'};
second_condition = {'B Cell'};
condition_names = [first_condition, second_condition];
first_data_threshold = 4;
second_data_threshold = 4;
dataset1_name_in_legend = 'H2';
dataset2_name_in_legend = 'M2';
series_matrix1_filename = 'newGSE107011_series_matrix';
series_matrix2_filename = 'newGSE124829_series_matrix';


exptable1_filename = '2human_pre_norm';
exptable2_filename = '2mouse_pre_norm';
title_ext = 'GSE107011(H2) and GSE124829(M2) before normalization';
filenames_ext = 'GSE107011_GSE124829_pre_norm';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)

exptable1_filename = '2human_cmaes_threshold';
exptable2_filename = '2mouse_cmaes_threshold';
title_ext = 'GSE107011(H2) and GSE124829(M2) after cmaes threshold=0.7';
filenames_ext = 'GSE107011_GSE124829_cmaes_threshold';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
                    

exptable1_filename = '2human_cmaes';
exptable2_filename = '2mouse_cmaes';
title_ext = 'GSE107011(H2) and GSE124829(M2) after cmaes with threshold=0';
filenames_ext = 'GSE107011_GSE124829_cmaes';
run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
                    
exptable1_filename = '2human_after_EB';
exptable2_filename = '2mouse_after_EB';
title_ext = 'GSE107011(H2) and GSE124829(M2) after EB';
filenames_ext = 'GSE107011_GSE124829_EB';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
                    
exptable1_filename = '2human_afterXPN';
exptable2_filename = '2mouse_afterXPN';
title_ext = 'between GSE107011(H2) and GSE124829(M2) after XPN';
filenames_ext = 'GSE107011_GSE124829_XPN';
run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
                    
% DEG_lists_from_paringNumber_to_gene_symbol('estimation_output\filtered_lists\',...
%                                     'GSE107011_GSE124829', 'translated_DEG_lists\')