clc

dir_name = 'estimation_output\filtered_lists\';

%H1_M1
disp('H1_M1 comparison')
calculate_diff_retaiment_idx(dir_name, 'H1_C1_H1_C2GSE60424_GSE122597_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'M1_C1_M1_C2GSE60424_GSE122597_pre_norm.csv')

%H1_M2
disp(newline)
disp('H1_M2 comparison')
calculate_diff_retaiment_idx(dir_name, 'H1_C1_H1_C2GSE60424_GSE124829_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'M2_C1_M2_C2GSE60424_GSE124829_pre_norm.csv')

%H2_M1
disp(newline)
disp('H2_M1 comparison')
calculate_diff_retaiment_idx(dir_name, 'H2_C1_H2_C2GSE107011_GSE122597_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'M1_C1_M1_C2GSE107011_GSE122597_pre_norm.csv')

%H2_M2
disp(newline)
disp('H2_M2 comparison')
calculate_diff_retaiment_idx(dir_name, 'H2_C1_H2_C2GSE107011_GSE124829_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'M2_C1_M2_C2GSE107011_GSE124829_pre_norm.csv')

%H1_H2
disp(newline)
disp('H1_H2 comparison')
calculate_diff_retaiment_idx(dir_name, 'H1_C1_H1_C2GSE60424_GSE107011_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'H2_C1_H2_C2GSE60424_GSE107011_pre_norm.csv')

%M1_M2
disp(newline)
disp('M1_M2 comparison')
calculate_diff_retaiment_idx(dir_name, 'M1_C1_M1_C2GSE122597_GSE124829_pre_norm.csv')
calculate_diff_retaiment_idx(dir_name, 'M2_C1_M2_C2GSE122597_GSE124829_pre_norm.csv')
