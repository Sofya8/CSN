function [samples_name, first_cell_type_logical_vector, sec_cell_type_logical_vector] ...
                    = filter_samples_and_clasify_groups(seriesTable,charName,valuesCharName)
    samples_name = filterByValues(seriesTable,seriesTable(1:end,1:1),charName,valuesCharName);
    first_cell_type_logical_vector = contains(samples_name(2,:),valuesCharName(1));
    sec_cell_type_logical_vector = contains(samples_name(2,:),valuesCharName(2));
end

function newCell = filterByValues(myCell,characteristics,var,values)
    myListOfValues = myCell(contains(characteristics,var),:); % extracts the values of the right characteristic
    neededsamples_idx = contains(myListOfValues, values);
    newCell = myCell(:,neededsamples_idx);
    newCell = [newCell(1, :);myListOfValues(neededsamples_idx)];
end