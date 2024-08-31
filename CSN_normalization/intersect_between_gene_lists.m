function [precent_genes_fst_table_in_intrsct,...
          precent_genes_sec_table_in_intrsct] = intersect_between_gene_lists(fst_table, sec_table)

    intrsct = intersect(fst_table ,sec_table);

    genes_num_fst_table = size(fst_table, 1);
    genes_num_sec_table = size(sec_table, 1);

    genes_num_intrsct = size(intrsct, 1);
    
    precent_genes_fst_table_in_intrsct = genes_num_intrsct / genes_num_fst_table;
    precent_genes_sec_table_in_intrsct = genes_num_intrsct / genes_num_sec_table;
end
