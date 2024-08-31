function calculate_diff_retaiment_idx(dir_name, filename)
    deg_list_name = strcat(dir_name, filename);
    pre_norm_gene_array = return_gene_array_of_lists(deg_list_name);
%     cmaes1_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'cmaes_human'));
%     cmaes2_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'cmaes_mouse'));
    baseline_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'baseline'));
    EB_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'EB'));
    DWD_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'DWD'));
    XPN_gene_array = return_gene_array_of_lists(strrep(deg_list_name, 'pre_norm', 'XPN'));

%     genes_arrays = {cmaes1_gene_array , cmaes2_gene_array };
    
%     method_names = {'cmaes1', 'cmaes2'};%
    method_names =  {'baseline','EB', 'DWD', 'XPN'};
    genes_arrays = {baseline_gene_array, EB_gene_array, DWD_gene_array, XPN_gene_array};
    generic_list_name = strrep(filename, '_pre_norm.csv', '');

    for i=1: 4
       C = intersect(pre_norm_gene_array, genes_arrays{i}); 
       intersection_size = size(C,1);
       C = union(pre_norm_gene_array, genes_arrays{i}); 
       union_size = size(C,1);
       diff_retaiment_idx = intersection_size/union_size;
       disp(strcat('The difference retaiment index for list   ', generic_list_name,...
                        ' with method ', method_names{i}, ' is: ', num2str(diff_retaiment_idx)))

    end
end

function gene_array = return_gene_array_of_lists(filename)
    table = readtable(filename);
    table = table(:,1);
    gene_array = table2array(table);
end