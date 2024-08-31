
function run_estimation(exptable1_filename, series_matrix1_filename, dataset1_name_in_legend,...
                        exptable2_filename, series_matrix2_filename, dataset2_name_in_legend,...
                        condition_field_name_in_series_matrix, first_condition, second_condition,...
                        first_data_threshold, second_data_threshold,...
                        filenames_ext, title_ext)
    
    build_folders_for_outputs()
    FDR_threshold = 0.05; 
    logFC_threshold = 1;
    
    fstExptable = readtable(strcat('data\', exptable1_filename));
    fstSeriesTable = readcell(strcat('data\', series_matrix1_filename));
    secExptable = readtable(strcat('data\', exptable2_filename));
    secSeriesTable = readcell(strcat('data\', series_matrix2_filename));
    
    [HC1_table, HC1_size] = create_group(fstExptable, fstSeriesTable,...
     condition_field_name_in_series_matrix, first_condition);
    [HC2_table, HC2_size] = create_group(fstExptable, fstSeriesTable,...
     condition_field_name_in_series_matrix, second_condition);
    [MC1_table, MC1_size] = create_group(secExptable, secSeriesTable,...
     condition_field_name_in_series_matrix, first_condition);
    [MC2_table, MC2_size] = create_group(secExptable, secSeriesTable,...
     condition_field_name_in_series_matrix, second_condition);
 
    [HC1_HC2, log_FC_HC1_HC2] = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               HC2_table, HC2_size, first_data_threshold, FDR_threshold, logFC_threshold);
	[MC1_MC2, log_FC_MC1_MC2]= find_DEGs(MC1_table, MC1_size, second_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold,logFC_threshold);
	[HC1_MC2, log_FC_HC1_MC2] = find_DEGs(HC1_table, HC1_size, first_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold);
	[MC1_HC2, log_FC_MC1_HC2] = find_DEGs(MC1_table, MC1_size, second_data_threshold,...
                               HC2_table, HC2_size, first_data_threshold, FDR_threshold, logFC_threshold);
    [HC1_MC1, log_FC_HC1_MC1] = find_DEGs(MC1_table, MC1_size, second_data_threshold,...
                                    HC1_table, HC1_size, first_data_threshold, FDR_threshold, logFC_threshold);
	[HC2_MC2, log_FC_HC2_MC2] = find_DEGs(HC2_table, HC2_size, first_data_threshold,...
                               MC2_table, MC2_size, second_data_threshold, FDR_threshold, logFC_threshold);
                           
	gene_lists = {HC1_HC2, MC1_MC2, HC1_MC2, MC1_HC2, HC1_MC1, HC2_MC2};
    log_FC_list = {log_FC_HC1_HC2, log_FC_MC1_MC2, log_FC_HC1_MC2, log_FC_MC1_HC2, log_FC_HC1_MC1,...
        log_FC_HC2_MC2};
    lists_filenames = {'H_C1_H_C2_list', 'M_C1_M_C2_list', 'H_C1_M_C2_list',...
                    'M_C1_H_C2_list', 'H_C1_M_C1_list','H_C2_M_C2_list'};
	intersects_figure_filename = strcat('intersects_heatmap');%, '.eps');
    intersects_figure_title = {...'Intersects between lists of differentially expressed genes', 
        title_ext};
    intersects_precent_figure_filename = strcat('intersects_heatmap_precent');%,'.eps');
    intersects_precent_figure_title = {...'Precent of differentially expressed genes from list A in intersects with list B', 
        title_ext};
    
    
    filter_and_illustrate_intersections_of_DEG(lists_filenames, gene_lists, intersects_figure_filename,...
        intersects_figure_title,intersects_precent_figure_filename, intersects_precent_figure_title,...
        title_ext, filenames_ext, dataset1_name_in_legend, dataset2_name_in_legend, log_FC_list);
end

function cellArray = readcell(path)
    T = readtable(path);
    cellArray = table2cell(T);
end

function [group_table, group_size] = create_group(exptable, seriesTable,...
    condition_field_name_in_series_matrix, condition)
    samples_name = filterByValues(seriesTable, seriesTable(:,1),...
                                condition_field_name_in_series_matrix, condition);   
    group_table = [exptable(:,1), exptable(:,samples_name(1,:))];
    group_size = size(group_table,2)-1;
    
end



