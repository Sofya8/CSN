function exptable_new = change_labels_of_genes(exptable_old, geneIDName_table, change_from_idx, to_idx )
    old_gene_labels = exptable_old(:,1);
    old_gene_labels.Properties.VariableNames(1) = {'PairingNumber'}; 
    lia = ismember (geneIDName_table(:,change_from_idx), old_gene_labels);
    geneIDName_table = geneIDName_table(lia,:);
    lia = ismember (old_gene_labels, geneIDName_table(:,change_from_idx));
    old_gene_labels = old_gene_labels(lia,:);
    exptable_old = exptable_old(lia,:);
    [~, locb] = ismember (geneIDName_table(:,change_from_idx), old_gene_labels);
    exptable_old = exptable_old(locb,:);
    new_lables = geneIDName_table(:,to_idx);
    exptable_temp = table2cell(exptable_old);
    new_lables = table2cell(new_lables);
    exptable_temp(:,1) = new_lables;
    exptable_new = cell2table(exptable_temp);  
    exptable_new.Properties.VariableNames(2:end) = exptable_old.Properties.VariableNames(2:end);
    exptable_new.Properties.VariableNames(1) = geneIDName_table.Properties.VariableNames(to_idx);
end