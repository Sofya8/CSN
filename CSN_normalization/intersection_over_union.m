function value = intersection_over_union(fst_table, sec_table)
    intrsct = intersect(fst_table ,sec_table);
    uni = union(fst_table ,sec_table);
    intrsct_size = size(intrsct, 1);
    uni_size = size(uni, 1);
    value = intrsct_size/uni_size;
end
    
    