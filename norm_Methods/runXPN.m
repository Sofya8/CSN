% Example: runXPN('ready_GSE124829','ready_GSE122597');

function runXPN(filename1,filename2)
    table1 = readtable(filename1);
    table2 = readtable(filename2);
    x1 = table2array(table1(:,2 :end));
    x2 = table2array(table2(:,2:end));
    [x1new, x2new] = XPN(x1,x2);
    x1new = array2table(x1new); 
    x2new = array2table(x2new);
    table1(:,2:end) = x1new;
    table2(:,2:end) = x2new;
    writetable(table1,strcat( strrep(filename1,'pre_norm',''), 'afterXPN'));
    writetable(table2,strcat( strrep(filename2,'pre_norm',''), 'afterXPN'));
end