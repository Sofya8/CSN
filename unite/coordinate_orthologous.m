addpath('../data/ortho_genes_tables')
addpath('../data/pre_processed_exp_tables')

fst_table_new_name = '../data/processed_exp_tables/1human_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/1mouse_pre_norm';
uniteByOne2OneOrthoGenes('ready_GSE60424', 'ready_GSE122597', fst_table_new_name, sec_table_new_name);

fst_table_new_name = '../data/processed_exp_tables/2human_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/2mouse_pre_norm';
uniteByOne2OneOrthoGenes('ready_GSE107011', 'ready_GSE124829', fst_table_new_name, sec_table_new_name);

fst_table_new_name = '../data/processed_exp_tables/3human_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/3mouse_pre_norm';
uniteByOne2OneOrthoGenes('ready_GSE60424', 'ready_GSE124829', fst_table_new_name, sec_table_new_name);

fst_table_new_name = '../data/processed_exp_tables/4human_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/4mouse_pre_norm';
uniteByOne2OneOrthoGenes('ready_GSE107011', 'ready_GSE122597', fst_table_new_name, sec_table_new_name);

fst_table_new_name = '../data/processed_exp_tables/mouse1_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/mouse2_pre_norm';
match_by_common_genes('ready_GSE122597', 'ready_GSE124829', fst_table_new_name, sec_table_new_name)

fst_table_new_name = '../data/processed_exp_tables/human1_pre_norm';
sec_table_new_name = '../data/processed_exp_tables/human2_pre_norm';
match_by_common_genes('ready_GSE60424', 'ready_GSE107011', fst_table_new_name, sec_table_new_name)
