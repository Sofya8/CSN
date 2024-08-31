function value = intersection_over_union(fst_table, sec_table)
    intrsct = intersect(fst_table ,sec_table);
    uni = union(fst_table ,sec_table);
    intrsct_size = size(intrsct);
    uni_size = size(uni);
    value = intrsct_size/uni_size;
end
    
    