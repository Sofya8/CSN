function DEG_lists_from_paringNumber_to_gene_symbol(input_files_dir, comparison_description, output_folder)
    
    output_dir = [output_folder comparison_description '\'];
    if ~exist(output_dir,'dir')
        mkdir(output_dir)
    end
    matfiles = dir(fullfile(strcat(input_files_dir, '*', comparison_description, '*.csv')));
    numberOfLists = size(matfiles,1);
    DEG_tables = cell(numberOfLists);
    DEG_table_names = cell(numberOfLists);
        for i = 1 : numberOfLists
           filename = matfiles(i).name;
           DEG_tables{i} = readtable([input_files_dir filename]);
           DEG_table_names{i} = filename;
        end
    for i=1:numberOfLists
        translatePairingNumbersInDEG_ListToGeneSymbol( DEG_tables{i},DEG_table_names{i}, output_dir);
    end
end
function translatePairingNumbersInDEG_ListToGeneSymbol(DEG_table, DEG_table_name, output_dir)
    human_label_table_name = 'new_ortho_human';
    mouse_label_table_name = 'new_ortho_mouse';
    human_label_table = readtable(human_label_table_name);
    mouse_label_table = readtable(mouse_label_table_name);
    DEG_table_with_names_mouse = change_labels_of_genes(DEG_table, mouse_label_table, 3, 2);
    DEG_table_with_names_human = change_labels_of_genes(DEG_table, human_label_table, 3, 2);
    DEG_table_with_names_mouse.Properties.VariableNames(1) = {'MouseGenes'};
    DEG_table_with_names_human.Properties.VariableNames(1) = {'HumanGenes'};
    new_DEG_table = [DEG_table_with_names_human(:,1) DEG_table_with_names_mouse];
    new_DEG_table_name = [output_dir DEG_table_name];
    writetable(new_DEG_table, new_DEG_table_name)
end