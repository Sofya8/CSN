function newCell = filterByValues(myCell,characteristics,var,values)
    myListOfValues = myCell(contains(characteristics,var),:); % extracts the values of the right characteristic
    neededsamples_idx = contains(myListOfValues, values);
    newCell = myCell(:,neededsamples_idx);
    newCell = [newCell(1, :);myListOfValues(neededsamples_idx)];
end