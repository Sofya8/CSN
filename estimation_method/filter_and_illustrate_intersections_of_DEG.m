function filter_and_illustrate_intersections_of_DEG(lists_file_names, list_of_DEG_tables, intersects_figure_filename,...
                                        intersects_figure_title, intersects_precent_figure_filename,...
                                        intersects_precent_figure_title,  title_ext, filenames_ext,...
                                        dataset1_name_in_legend, dataset2_name_in_legend,...
                                        log_FC_list)
    
	lists_file_names = change_lists_file_names(lists_file_names, dataset1_name_in_legend, dataset2_name_in_legend);
                                                
    save_tables(list_of_DEG_tables, lists_file_names, 'estimation_output\filtered_lists\', filenames_ext);
    
    [list_of_separated_tables,...
        separated_tables_names] = separate_lists_by_log_fold_change_values(list_of_DEG_tables,...
                                                                              lists_file_names, log_FC_list);
    
    [intersectMatrix, proportionMatrix] = intersects_between_lists(list_of_DEG_tables,...
                            lists_file_names, 'estimation_output\intersections\', filenames_ext);
                        
	value = calculate_quantative_value(proportionMatrix);
    display(strcat("The CSC value between ", title_ext, " is: ", num2str(value)));
    
    
    [separated_intersectMatrix_pre_norm, separated_precentMatrix_pre_norm] = ...
                        intersects_between_lists(list_of_separated_tables,...
                        separated_tables_names, 'estimation_output\separated_intersections\',...
                        filenames_ext);
    separated_legend_names = strrep(separated_tables_names,'_lower', '<');
    separated_legend_names = strrep(separated_legend_names,'_higher', '>');
 
    fignum = 1;
    fignum = create_and_save_heatmap_for_matrix(intersectMatrix, lists_file_names, intersects_figure_title,...
                                                        strcat('estimation_output\', ...
                                                                intersects_figure_filename,...
                                                                filenames_ext),...
                                                        fignum);
	fignum = create_and_save_heatmap_for_matrix(separated_intersectMatrix_pre_norm, separated_legend_names, intersects_figure_title,...
                                                        strcat('estimation_output\', 'separated_lists',...
                                                                intersects_figure_filename,...
                                                                filenames_ext), fignum);
% 	fignum = create_and_save_heatmap_for_matrix(proportionMatrix, lists_file_names, intersects_precent_figure_title,...
%                                                         strcat('estimation_output\',...
%                                                                 intersects_precent_figure_filename,...
%                                                                 filenames_ext),...
%                                                         fignum);                                          
%     create_and_save_heatmap_for_matrix(separated_precentMatrix_pre_norm, separated_legend_names, intersects_precent_figure_title,...
%                                                        strcat('estimation_output\', 'separated_lists',...
%                                                                 intersects_precent_figure_filename,...
%                                                                 filenames_ext),...
%                                                        fignum);
end

function legend_names = change_lists_file_names(lists_file_names, dataset1_name_in_legend, dataset2_name_in_legend)
    
    lists_file_names = strrep(lists_file_names,'_list','');
    lists_file_names = strrep(lists_file_names, 'H_', strcat(dataset1_name_in_legend,'_'));
    legend_names = strrep(lists_file_names, 'M_', strcat(dataset2_name_in_legend,'_'));
end

function value = calculate_quantative_value(precentMatrix)
    value1 = (sum(sum(precentMatrix(1:4,1:4)))-4)/12;
    value2 = (sum(sum(precentMatrix(1:4,5:6))) + sum(sum(precentMatrix(5:6,1:4))))/16;
    value = value1/value2;
end

function [list_of_separated_tables,...
            separated_tables_names] = separate_lists_by_log_fold_change_values(list_of_tables,...
                                                                        tables_names, log_FC_list)
    number_of_lists = length(tables_names);
    list_of_tables_fst_bigger = cell(1,number_of_lists);
    list_of_tables_sec_bigger = cell(1,number_of_lists);
    for i = 1:number_of_lists
        table_i = list_of_tables{i};
        list_of_tables_fst_bigger{i} = table_i(log_FC_list{i}<0,:);
        list_of_tables_sec_bigger{i} = table_i(log_FC_list{i}>0,:);
    end
    
    list_of_separated_tables = [list_of_tables_fst_bigger, list_of_tables_sec_bigger];
    separated_tables_names = [strcat(tables_names, '_higher'), strcat(tables_names, '_lower')];
end

function [intersectMatrix, precentMatrix] = intersects_between_lists(list_of_tables, tables_names,...
                                                                        dir_name, filenames_ext)
    number_of_lists = length(tables_names);
    intersectMatrix = zeros(number_of_lists,number_of_lists);
    precentMatrix = zeros(number_of_lists,number_of_lists);
    for i = 1 : number_of_lists
       table_i = list_of_tables{i};
       intersectMatrix(i, i) = size(table_i,1);
       precentMatrix(i, i) = 1;
       for j = i+1 : number_of_lists
          [genes_num_intrsct,precent_genes_fst_table_in_intrsct,...
                  precent_genes_sec_table_in_intrsct] = intersect_between_gene_lists(list_of_tables{i},...
                                                            list_of_tables{j}, tables_names{i},...
                                                            tables_names{j}, dir_name, filenames_ext);
          intersectMatrix(i, j) = genes_num_intrsct;
          intersectMatrix(j, i) = genes_num_intrsct;
          precentMatrix(i, j) = precent_genes_fst_table_in_intrsct;
          precentMatrix(j, i) = precent_genes_sec_table_in_intrsct;
       end
    end
    precentMatrix(isnan(precentMatrix))=0;
end

function list_of_tables = build_array_of_filtered_tables(tables_names, FDR_threshold, logFC_threshold)

    number_of_lists = length(tables_names);
    
    list_of_tables = cell(1,number_of_lists);
    for i = 1:number_of_lists
        list_of_tables{i} = filter_differentially_expressed_genes(tables_names{i}, FDR_threshold, logFC_threshold); 
    end
end

function save_tables(list_of_tables, DEG_table_name, dir, filenames_ext)
    
    number_of_lists = length(list_of_tables);
    for i = 1:number_of_lists
       writetable(list_of_tables{i} ,strcat(dir, DEG_table_name{i}, filenames_ext, '.csv')); 
    end 
end

function [genes_num_intrsct, precent_genes_fst_table_in_intrsct,...
    precent_genes_sec_table_in_intrsct] = intersect_between_gene_lists(fst_table, sec_table,...
                                            fst_table_name, sec_table_name, dir_name, filenames_ext)

    intrsct = intersect(fst_table.tracking_id ,sec_table.tracking_id);
    intrsct = array2table(intrsct);
    if isempty(intrsct)
        intrsct = table('Size', [0, 1], 'VariableTypes', {'double'}, 'VariableNames', {'tracking_id'});
    else
        intrsct.Properties.VariableNames = {'tracking_id'};
    end
    writetable(intrsct ,strcat(dir_name, 'intersect_', fst_table_name,...
                                            '_and_', sec_table_name, filenames_ext, '.csv'));

    genes_num_fst_table = size(fst_table, 1);
    genes_num_sec_table = size(sec_table, 1);

    genes_num_intrsct = size(intrsct, 1);
    
    precent_genes_fst_table_in_intrsct = genes_num_intrsct / genes_num_fst_table;
    precent_genes_sec_table_in_intrsct = genes_num_intrsct / genes_num_sec_table;
end

function table = filter_differentially_expressed_genes(DEG_table_name, FDR_threshold, logFC_threshold)
    table = readtable(DEG_table_name);
    table(table.FDR > FDR_threshold ,:) = [];
    table(abs(table.logFC) < logFC_threshold ,:) = [];
end

function fignum = create_and_save_heatmap_for_matrix(matrix, list_of_names, plot_title,...
                                                        plot_file_name, fignum)
    h = figure;
    list_of_names = strrep(list_of_names,'_C','C');
    list_of_names = strrep(list_of_names,'_C','C');
    create_heatmap(list_of_names,matrix,plot_title, h)
    %saveas(h ,plot_file_name, 'jpg'); %%for saving the plots
    fignum = fignum + 1;
end
function create_heatmap(list_of_names,matrix,plot_title, h1)
     h = heatmap(h1, list_of_names,list_of_names,matrix);
     h.ColorbarVisible = 'off';
     if ~strcmp(plot_title,'')
        title(plot_title);
     end
end