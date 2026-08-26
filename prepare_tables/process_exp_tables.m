fig = figure;
addpath('../data/raw_exp_tables')
addpath('../data/series')
addpath('../data/raw_series')
addpath('../data/pre_processed_exp_tables')
processTableGSE60424()
processTableGSE107011()
processTableGSE122597()
processTableGSE124829()




function processTableGSE107011()
    expTablePath = 'GSE107011_Processed_data_TPM.txt';
    exptable = readtable(expTablePath);
    
    genes = table2cell(exptable(:,1));
    with_dot = contains(genes,'.');
    genes(with_dot)= extractBefore(genes(with_dot),".");
    exptable(:,1) = cell2table(genes);
    namesTablePath = 'GSE107011_names';
    names_table = readtable(namesTablePath);

    names_table.Sample_title = extractBefore(names_table.Sample_title, "_rep");
    names_table.Sample_title = strrep(names_table.Sample_title, '9', 'x9');
    exptable = rename_by_names_file(exptable, names_table);
    exptable.Properties.VariableNames(1) = {'tracking_id'};
    filtered_table = exptable(:,~contains(exptable.Properties.VariableNames,'Neutrophils'));
    normalized_exptable = exptable_to_log_normalized_scale(filtered_table);
    writetable(normalized_exptable,'../data/pre_processed_exp_tables/ready_GSE107011');
    plot_box_plots(exptable, normalized_exptable, 2, '(GSE107011 - H2)');
end

function processTableGSE122597()
    expTablePath = 'GSE122597_Gene_count_table.csv';
    exptable = readtable(expTablePath);
    exptableNames = exptable.Properties.VariableNames;
    geneID_name_table = readtable('mouse_geneID_name');
    exptable = table2cell(exptable);
    geneID_name_table = table2cell(geneID_name_table);
    exptable = change_labels_of_genes(exptable, geneID_name_table , 2, 1);
    exptable = cell2table(exptable);
     exptable.Properties.VariableNames = ['tracking_id', exptableNames(2:end)];
    normalized_exptable = exptable_to_log_normalized_scale(exptable);
    writetable(normalized_exptable,'../data/pre_processed_exp_tables/ready_GSE122597');
    plot_box_plots(exptable, normalized_exptable, 3, '(GSE122597 - M1)')
end


function processTableGSE60424()
    expTablePath = 'GSE60424_GEOSubmit_FC1to11_normalized_counts.txt';
    exptable = readtable(expTablePath);
    
    
    genes = table2cell(exptable(:,1));
    with_dot = contains(genes,'.');
    genes(with_dot)= extractBefore(genes(with_dot),".");
    exptable(:,1) = cell2table(genes);

    namesTablePath = 'GSE60424_names';
    names_table = readtable(namesTablePath);
    names_table.name = strrep(names_table.name, '-', '_'); 

    exptable = rename_by_names_file(exptable, names_table); 
    keep = [true, contains(exptable.Properties.VariableNames(2:end), 'Healthy_Control')];
    exptable = exptable(:, keep);
    exptable = exptable(:,~contains(exptable.Properties.VariableNames,'Neutrophils' )& ~contains(exptable.Properties.VariableNames, 'Whole'));
    exptable.Properties.VariableNames(1) = {'tracking_id'};
    normalized_exptable = exptable_to_log_normalized_scale(exptable);
    plot_box_plots(exptable, normalized_exptable, 1, '(GSE60424 - H1)')
    writetable(normalized_exptable,'../data/pre_processed_exp_tables/ready_GSE60424');
end


function processTableGSE124829()
    filtered_seriestable = readcell('newGSE124829_series_matrix');
    expTablePath = 'GSE124829_Gene_count_table.csv';
    exptable = readtable(expTablePath);
    geneID_name_table = readtable('mouse_geneID_name');
    sample_names = exptable.Properties.VariableNames;
    exptable_without_sample_names = table2cell(exptable);
    geneID_name_table = table2cell(geneID_name_table);
    exptable_without_sample_names = change_labels_of_genes(exptable_without_sample_names, geneID_name_table , 2, 1);
    exptable = cell2table(exptable_without_sample_names);
    exptable.Properties.VariableNames = sample_names;
    exptable.Properties.VariableNames(1) = {'tracking_id'};
    filterd_exptableNames = filtered_seriestable(1,2:end);
    exptable = [exptable(:,1), exptable(:,filterd_exptableNames)];
    exptable = [exptable(:,1:13), exptable(:,15:end)];
    normalized_exptable = exptable_to_log_normalized_scale(exptable);
    writetable(normalized_exptable,'../data/pre_processed_exp_tables/ready_GSE124829');
    plot_box_plots(exptable, normalized_exptable, 4, '(GSE124829 - M2)')
end

function exptable = exptable_to_log_normalized_scale(exptable)
	counts = table2array(exptable(:,2:end));
    normCounts = normalize_read_counts_to_library_size(counts);
    normCounts(normCounts<1) = 1;
    log_normCounts = log2(normCounts);
    exptable(:,2:end) = array2table(log_normCounts);
end

function normCounts = normalize_read_counts_to_library_size(counts)
    pseudoRefSample = geomean(counts,2);
    nz = pseudoRefSample > 0;
    ratios = bsxfun(@rdivide,counts(nz,:),pseudoRefSample(nz));
    sizeFactors = median(ratios,1);
    normCounts = bsxfun(@rdivide,counts,sizeFactors);
end

function exptable = rename_by_names_file(exptable, names_table)

    % Get the expression-table sample columns (everything except genenames)
    sampleCols = exptable.Properties.VariableNames(2:end);

    % Match each exptable column to Sample_title
    [tf, idx] = ismember(sampleCols, names_table.Sample_title);

    % Keep only columns that were found in names_table
    exptable = exptable(:, [true, tf]);

    % Rename the matched columns using names_table.name
    newNames = names_table.name(idx(tf));

    exptable.Properties.VariableNames(2:end) = newNames;

    
end 

function plot_box_plots(counts, normCounts, exp_num, exp_name)
    counts = table2array(counts(:,2:end));
    normCounts = table2array(normCounts(:,2:end));


    subplot(2,2,exp_num)
    bins = 0:1:15;
    h = hist(normCounts, bins);
    plot(bins, h);
    title(strcat('Normalized Counts ', exp_name));
    xlim([1 15])
    
%     subplot(2,2,exp_num)
%     bins = 0:1:15;
%     h = hist(log2(counts), bins);
%     plot(bins, h);
%     title(strcat('Non-normalized Counts ', exp_name));
%     xlim([1 15])
    
    
%     subplot(2,4,exp_num)
%     maboxplot(log2(counts),'title',strcat('Counts ', exp_name),'orientation','horizontal')
%     ylabel('sample')
%     xlabel('log2(counts)')
% 
%     subplot(2,4,4+exp_num)
%     maboxplot(normCounts,'title',strcat('Normalized Counts ', exp_name),'orientation','horizontal')
%     ylabel('sample')
%     xlabel('log2(normCounts)')
end