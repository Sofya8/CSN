
function log_fold_change_array = calc_log_FC(fst_group_array, sec_group_array)
    fst_group_mean = mean(fst_group_array, 2);
    sec_group_mean = mean(sec_group_array, 2);
    log_fold_change_array = fst_group_mean - sec_group_mean;
    
end