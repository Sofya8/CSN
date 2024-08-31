function [filtered_fst_group_table, filtered_sec_group_table] = filter_genes_under_threshold(...
                        fst_group_table, fst_group_size, fst_dataset_threshold,...
                        sec_group_table, sec_group_size, sec_dataset_threshold)
    both_groups_array = table2array([fst_group_table(:,2:end), sec_group_table(:,2:end)]);
    fst_threshold_vec = fst_dataset_threshold*ones(1, fst_group_size);
    sec_threshold_vec = sec_dataset_threshold*ones(1, sec_group_size);
    fst_idx_matrix = both_groups_array(:,1:fst_group_size) > fst_threshold_vec;
    sec_idx_matrix = both_groups_array(:,1+fst_group_size:end) > sec_threshold_vec;
    overall_idx = sum([fst_idx_matrix sec_idx_matrix], 2)<min(fst_group_size,sec_group_size);
    filtered_fst_group_table = fst_group_table(~overall_idx, :);
    filtered_sec_group_table = sec_group_table(~overall_idx, :);
    
end