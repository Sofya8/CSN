function baseline(filename1,filename2)
    table1 = readtable(filename1);
    table2 = readtable(filename2);
    x1 = table2array(table1(:,2:end));
    x2 = table2array(table2(:,2:end));
    [x1new, x2new]  = GetMatrixByOrderOfTheMean(x1,x2);
    x1new = array2table(x1new); 
    x2new = array2table(x2new);
    table1(:,2:end) = x1new;
    table2(:,2:end) = x2new;
    writetable(table1,strcat(strrep(filename1,'pre_norm',''), 'after_baseline'));
    writetable(table2,strcat(strrep(filename2,'pre_norm',''), 'after_baseline'));
end


function [new_matrix1, new_matrix2] = GetMatrixByOrderOfTheMean(matrix1, matrix2)
    %% Dor Marciano's code
    values = mean(matrix1,2);
    values = sort(values);
    matrix = [matrix1 matrix2];
    matrix1_sample_size = size(matrix1, 2);
    valuesMatrix = repmat(values, 1, size(matrix,2)); 
    [~, sortIndices] = sort(matrix); 
    whichColumnMatrix = repmat(1:size(matrix,2), size(matrix,1), 1);
    vectorOfLinearIndicesInCorrectOrder = sub2ind(size(matrix), sortIndices(:), whichColumnMatrix(:));
    valuesMatrixCopy = valuesMatrix;
    valuesMatrixCopy(vectorOfLinearIndicesInCorrectOrder) = valuesMatrix(:);
    new_matrix1 = valuesMatrixCopy(:,1:matrix1_sample_size);
    new_matrix2 = valuesMatrixCopy(:,1+matrix1_sample_size:end);
end