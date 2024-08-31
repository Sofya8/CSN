function match_by_common_genes(filename1, filename2, new_filename1, new_filename2)
    table1 = readtable(filename1);
    table2 = readtable(filename2);
    [~, ia, ib] = intersect(table1(:,1), table2(:,1));
    table1 = table1(ia,:);
    table2 = table2(ib,:);
    writetable(table1,new_filename1);
    writetable(table2,new_filename2);
end