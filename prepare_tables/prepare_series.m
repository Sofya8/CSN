addpath('../data/raw_series')

processSeriesTableGSE107011();
processSeriesTableGSE122597();
processSeriesTableGSE60424();
processSeriesTableGSE124829()

function processSeriesTableGSE107011()
%% TODO: take care of read and write cell
%% TODO: take care of real series matrix(without deleting first 20 lines by hand)
    filename = 'GSE107011_series_matrix';
    myCell = readcell(filename);
    myCell(10,1) = {'!Sample_characteristics_name'};
    %myCell = [myCell(:,1) filterByValue(myCell,myCell(:,1),'!Sample_characteristics_name','Naive')];
    myCell = unify_B_cells_names(myCell);
    myCell(1,:) = strrep(myCell(19,:), '9', 'x9');
    
    namesTablePath = 'GSE107011_names';
    names_table = readcell(namesTablePath);
    names_table(:,3) = extractBefore(names_table(:,3), "_rep"); 
    names_table(:,3) = strrep(names_table(:,3),  '9', 'x9'); 
    
    rotated_cell_matrix = myCell(:,2:end)';
    rotated_cell_matrix = change_labels_by_lists(rotated_cell_matrix, names_table, 3, 2);
    myCell = [myCell(:,1) rotated_cell_matrix'];
    
    writecellarray(myCell, '../data/series/newGSE107011_series_matrix');
end

function processSeriesTableGSE122597()
    filename = 'GSE122597_series_matrix';
    myCell = readcell(filename);
    myCell(20,1) = {'!Sample_characteristics_name'};
    myCell = unify_B_cells_names(myCell);
    myCell(1,:) = strrep(myCell(1,:), '.', '_');
    myCell(1,:) = strrep(myCell(1,:), '#', '_');
    writecellarray(myCell, '../data/series/newGSE122597_series_matrix');
end


function processSeriesTableGSE124829()
%% TODO: take care of read and write cell
%% TODO: take care of real series matrix(without deleting first 20 lines by hand)
    filename = 'GSE124829_series_matrix';
   
    dataset_partition_table = readtable('GSE124829_partition.xlsx');
   % exptable = exptable(B_dataset_sample_names,:);
    
    myCell = readcell(filename);
    myCell(22,1) = {'!Sample_characteristics_name'};
    B_dataset_entries = contains(dataset_partition_table.Dataset, 'B');
   % B_dataset_sample_names = dataset_partition_table.x_Sample_geo_accession(B_dataset_entries);
    samples_only = myCell(:,2:end);
    samples_only = samples_only(:,B_dataset_entries);
    myCell =  [myCell(:,1) samples_only];
    treatment_samples = contains(myCell(22,:), 'IFN');
    myCell(:, treatment_samples) = [];
    myCell = unify_B_cells_names(myCell);
    myCell(1,:) = strrep(myCell(1,:), '.', '_');
    myCell(1,:) = strrep(myCell(1,:), '#', '_');
    writecellarray(myCell, '../data/series/newGSE124829_series_matrix');
end


function processSeriesTableGSE60424()
%% TODO: take care of read and write cell
%% TODO: take care of real series matrix(without deleting first 20 lines by hand)
    filename = 'GSE60424_series_matrix';
    myCell = readcell(filename);
    myCell(12,1) = {'!Sample_characteristics_name'};
     myCell(14,1) = {'!Sample_diseasestatus_ch1'};
     myCell = [myCell(:,1) filterByValue(myCell,myCell(:,1),'!Sample_diseasestatus_ch1','diseasestatus: Healthy Control')];
    myCell = unify_B_cells_names(myCell);
    
    namesTablePath = 'GSE60424_names';
    names_table = readcell(namesTablePath);
    names_table(:,2) = strrep(names_table(:,2), '-', '_'); 

    rotated_cell_matrix = myCell(:,2:end)';
    rotated_cell_matrix = change_labels_by_lists(rotated_cell_matrix, names_table, 4, 2);
    myCell = [myCell(:,1) rotated_cell_matrix'];
    myCell(:,contains(myCell(1,:), 'Whole_Blood')) = [];

    
    writecellarray(myCell, '../data/series/newGSE60424_series_matrix');
end

function cellArray = readcell(path)
    T = readtable(path, 'Delimiter', '\t');
    cellArray = table2cell(T);
end

function writecellarray(cellarray, path)
    table = cell2table(cellarray);
    table.Properties.VariableNames = [{'Sample_title'} table2cell(table(1,2:end))];
    table(1:end-1,:) = table(2:end,:);
    table(end,:)=[];
    writetable(table, path);
end
function newCell = filterByValue(myCell,characteristics,var,value)
    myListOfValues = myCell(contains(characteristics,var),:); % extracts the values of the right characteristic
    newCell = myCell(:,contains(myListOfValues, value)); % filters the samples of the value 
end

function myCell = unify_B_cells_names(myCell)
    myCell = strrep(myCell, 'B-cells', 'B Cell');
    myCell = strrep(myCell, 'B cells', 'B Cell');
end