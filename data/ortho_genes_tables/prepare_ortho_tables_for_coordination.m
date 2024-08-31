function prepare_ortho_tables_for_coordination

%% This is how the ortho_gene tables from ensemble were processed

    tablehuman = readtable('original_ortho-homo');
    tablemice = readtable('original_ortho-mouse');
    num_of_genes = size(tablehuman, 1);
    col = linspace(1,num_of_genes,num_of_genes)';
    tablehuman.PairingNumber = col;
    tablemice.PairingNumber = col;
    writetable(tablehuman, 'new_ortho_human');
    writetable(tablemice, 'new_ortho_mouse');
end