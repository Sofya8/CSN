function [group_table, group_size] = create_group(exptable, seriesTable,...
    condition_field_name_in_series_matrix, condition)
    samples_name = filterByValues(seriesTable, seriesTable(:,1),...
                                condition_field_name_in_series_matrix, condition);   
    group_table = [exptable(:,1), exptable(:,samples_name(1,:))];
    group_size = size(group_table,2)-1;
    
end