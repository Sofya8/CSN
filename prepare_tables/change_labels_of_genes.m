function exptable = change_labels_of_genes(exptable, geneIDName_table, change_from_idx, to_idx )
    old_gene_labels = exptable(:,1);
    lia = ismember (geneIDName_table(:,change_from_idx), old_gene_labels);
    geneIDName_table = geneIDName_table(lia,:);
    lia = ismember (old_gene_labels, geneIDName_table(:,change_from_idx));
    old_gene_labels = old_gene_labels(lia,:);
    exptable = exptable(lia,:);
    [~, locb] = ismember (geneIDName_table(:,change_from_idx), old_gene_labels);
    exptable = exptable(locb,:);
    exptable(:,1) = geneIDName_table(:,to_idx);
end