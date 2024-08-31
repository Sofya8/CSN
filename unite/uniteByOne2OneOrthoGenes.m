function table = uniteByOne2OneOrthoGenes(human_exptable_filename, mouse_exptable_filename,...
                                          human_table_new_filename, mouse_exptable_new_filename)
    mouse_exp_table = readtable(mouse_exptable_filename);
    human_exp_table = readtable(human_exptable_filename);
    tablehuman = readtable('new_ortho_human');
    tablemice = readtable('new_ortho_mouse');
    [tablehuman,tablemice] = extract_one2one_orthologous(tablehuman,tablemice);
    [tablehuman,tablemice] = existedOrthoGenesForExperiment(tablehuman,tablemice,human_exp_table);
    [tablemice, tablehuman] = existedOrthoGenesForExperiment(tablemice,tablehuman,mouse_exp_table);
    tablehuman_gene_id = tablehuman(:,1);
    tablemice_gene_id = tablemice(:, 1);
    lia = ismember(mouse_exp_table(:,1),tablemice_gene_id);
    mouse_exp_table = mouse_exp_table(lia,:);
    [~,locb] = ismember(tablemice_gene_id, mouse_exp_table(:,1));
    mouse_exp_table = mouse_exp_table(locb,:);
    lia = ismember(human_exp_table(:,1),tablehuman_gene_id);
    human_exp_table = human_exp_table(lia,:);
    [~,locb] = ismember(tablehuman_gene_id, human_exp_table(:,1));
    human_exp_table = human_exp_table(locb,:);
    table = [mouse_exp_table, human_exp_table(:,2:end)]; 
    
    
    mouse_exp_table.tracking_id = tablemice.PairingNumber;
    human_exp_table.tracking_id = tablehuman.PairingNumber;
    
    writetable(mouse_exp_table, mouse_exptable_new_filename);
    writetable(human_exp_table, human_table_new_filename);
end

function [fst_species_table, sec_species_table] = extract_one2one_orthologous(fst_species_table, sec_species_table)
    fst_species_reapetedRows = findRepeatedRowsOfTable(fst_species_table);
    sec_species_reapetedRows = findRepeatedRowsOfTable(sec_species_table);
    sum_of_repeated_values_indecses = fst_species_reapetedRows + sec_species_reapetedRows;
    reapetedRows = sum_of_repeated_values_indecses > 0;
    fst_species_table(reapetedRows,:) = [];
    sec_species_table(reapetedRows, :) = [];
end

function reapetedRows = findRepeatedRowsOfTable(geneID_table)
    geneID_array = table2cell(geneID_table(:,1));
    [~,gene_indeces_without_copies] = unique(geneID_array);
    ortho_many_genes = geneID_array;
    ortho_many_genes(gene_indeces_without_copies) = [];
    reapetedRows = contains(geneID_array, ortho_many_genes);
end

function [tableFstSpecie,tableSecSpecie] = existedOrthoGenesForExperiment(tableFstSpecie,tableSecSpecie,tableExperiment) 
    lia = ismember(tableFstSpecie(:,1), tableExperiment(:,1));
    tableFstSpecie = tableFstSpecie(lia,:);
    tableSecSpecie = tableSecSpecie(lia,:);
end