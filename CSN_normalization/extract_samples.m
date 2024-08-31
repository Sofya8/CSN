% function [samples_name, fst_group_size, sec_group_size]  = extract_samples(seriesTable,condition_field_name_in_series_matrix,...
%                                                             condition_names)

function [fst_samples_name, sec_samples_name]  = extract_samples(seriesTable,condition_field_name_in_series_matrix,...
                                                            condition_names)

    [samples_name, first_cell_type_logical_vector, sec_cell_type_logical_vector] = filter_samples_and_clasify_groups(seriesTable,condition_field_name_in_series_matrix,condition_names);
    fst_samples_name = samples_name(:,first_cell_type_logical_vector);
    sec_samples_name = samples_name(:,sec_cell_type_logical_vector);
%     samples_name = [fst_samples_name sec_samples_name];
%     fst_group_size = sum(first_cell_type_logical_vector);
%     sec_group_size = sum(sec_cell_type_logical_vector);
    
end
