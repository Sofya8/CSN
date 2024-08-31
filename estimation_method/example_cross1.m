clc
close all

condition_field_name_in_series_matrix = '!Sample_characteristics_name';
first_condition = {'CD4'};
second_condition = {'B Cell'};
condition_names = [first_condition, second_condition];
first_data_threshold = 4;
second_data_threshold = 4;
dataset1_name_in_legend = 'H1';
dataset2_name_in_legend = 'M1';
series_matrix1_filename = 'newGSE60424_series_matrix';
series_matrix2_filename = 'newGSE122597_series_matrix';



exptable1_filename = '1human_pre_norm';
exptable2_filename = '1mouse_pre_norm';
title_ext = 'GSE60424(H1) and GSE122597(M1) before normalization';
filenames_ext = 'GSE60424_GSE122597_pre_norm';
run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)


exptable1_filename = '1human_cmaes_threshold';
exptable2_filename = '1mouse_cmaes_threshold';
title_ext = 'GSE60424(H1) and GSE122597(M1) after CMAES with threshold=0.7';
filenames_ext = 'GSE60424_GSE122597_cmaes_threshold';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)

exptable1_filename = '1human_cmaes';
exptable2_filename = '1mouse_cmaes';
title_ext = 'GSE60424(H1) and GSE122597(M1) after CMAES with threshold=0';
filenames_ext = 'GSE60424_GSE122597_cmaes';


run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)

exptable1_filename = '1human_after_EB';
exptable2_filename = '1mouse_after_EB';
title_ext = 'GSE60424(H1) and GSE122597(M1) after EB';
filenames_ext = 'GSE60424_GSE122597_EB';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)

exptable1_filename = '1human_afterXPN';
exptable2_filename = '1mouse_afterXPN';
title_ext = 'GSE60424(H1) and GSE122597(M1) after XPN';
filenames_ext = 'GSE60424_GSE122597_XPN';

run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename,dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
                    
% DEG_lists_from_paringNumber_to_gene_symbol('estimation_output\filtered_lists\',...
%                                     'GSE60424_GSE122597', 'translated_DEG_lists\')
