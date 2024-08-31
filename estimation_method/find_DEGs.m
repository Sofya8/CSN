function [gene_list, log_FC_list] = find_DEGs(fst_group_table, fst_group_size, fst_dataset_threshold,...
                               sec_group_table, sec_group_size, sec_dataset_threshold, FDR_threshold,...
                               log_fold_change_threshold)
    
    [filtered_fst_group_table, filtered_sec_group_table] = filter_genes_under_threshold(...
                        fst_group_table, fst_group_size, fst_dataset_threshold,...
                        sec_group_table, sec_group_size, sec_dataset_threshold);

    fst_group_array = table2array(filtered_fst_group_table(:, 2:end));
    sec_group_array = table2array(filtered_sec_group_table(:, 2:end));
    %[h,PValues,ci,stats] = ttest2(fst_group_array', sec_group_array', 'Vartype', 'unequal');
    [PValues, TScores] = mattest(fst_group_array, sec_group_array);%,'showhist',true,'showplot',true);
    FDR = mafdr(PValues, 'BHFDR',true);

%     tLocal = nbintest(fst_group_array,sec_group_array, 'VarianceLink','LocalRegression');
%     FDR = mafdr(tLocal.pValue,'BHFDR',true);
    
    indecies_of_significant_genes_by_t_test = FDR<FDR_threshold;
    
    log_fold_change_array = calc_log_FC(fst_group_array, sec_group_array);
    indecies_by_log_fold_change = abs(log_fold_change_array) > log_fold_change_threshold;
    
    indecies_of_significant_genes = (indecies_by_log_fold_change +...
                                       indecies_of_significant_genes_by_t_test)==2;
    
    gene_list = filtered_fst_group_table(indecies_of_significant_genes,1);
    log_FC_list = log_fold_change_array(indecies_of_significant_genes);
    
end