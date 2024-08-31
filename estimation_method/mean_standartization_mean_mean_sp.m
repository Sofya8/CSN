condition_field_name_in_series_matrix = '!Sample_characteristics_name';
condition = {'CD4'};


H1_series_matrix_filename = strcat('data\','newGSE60424_series_matrix');
H1_seriesTable = readcell(H1_series_matrix_filename);
% H2_series_matrix_filename = strcat('data\','newGSE107011_series_matrix');
% H2_seriesTable = readcell(H2_series_matrix_filename);
M1_series_matrix_filename = strcat('data\','newGSE122597_series_matrix');
M1_seriesTable = readcell(M1_series_matrix_filename);
% M2_series_matrix_filename = strcat('data\','newGSE124829_series_matrix');
% M2_seriesTable = readcell(M2_series_matrix_filename);

X_axis_series = H1_seriesTable;%, {H2_seriesTable}, {H1_seriesTable}...
                  % {H2_seriesTable}, {H1_seriesTable}, {M1_seriesTable}];
               
Y_axis_series = M1_seriesTable;%, {M2_seriesTable}, {M2_seriesTable}...
                   %{M1_seriesTable}, {H2_seriesTable}, {M2_seriesTable}];
X_axis_exptable_name = '1human_pre_norm';%, {'2human_pre_norm'}, {'3human_pre_norm'},...
                     %{'4human_pre_norm'}, {'human1_pre_norm'}, {'mouse1_pre_norm'}];
                     
X_axis_exptable_name = strcat('data\',X_axis_exptable_name);
                        
                        

Y_axis_exptable_name = '1mouse_pre_norm';%, {'2mouse_pre_norm'}, {'3mouse_pre_norm'},...
%                      {'4mouse_pre_norm'}, {'human2_pre_norm'}, {'mouse2_pre_norm'}];
                     
Y_axis_exptable_name = strcat('data\',Y_axis_exptable_name);


origin_tableX = readtable(X_axis_exptable_name);
origin_tableY = readtable(Y_axis_exptable_name);
stdized_table_X = origin_tableX;
stdized_table_Y = origin_tableY;

temp_X = table2array(origin_tableX(:,2:end));
temp_Y = table2array(origin_tableY(:,2:end));

median_X = median(temp_X,2);
median_Y = median(temp_Y,2);
temp_X = temp_X + 0.5*(median_Y - median_X);
temp_Y = temp_Y + 0.5*(median_X - median_Y);

temp_X = array2table(temp_X);
temp_Y = array2table(temp_Y);
stdized_table_X(:,2:end) = temp_X(:,:);
stdized_table_Y(:,2:end) = temp_Y(:,:);
% 
% exptables = [{origin_tableX },{origin_tableY}, {stdized_table_X}, stdized_table_Y}];
% seriesTables = [{H1_seriesTable}, {M1_seriesTable}]

vector_X_axis = {};
vector_Y_axis = {};
vector_X_axis{1} = extract_mean_vector_of_condition(origin_tableX, H1_seriesTable,...
                                                condition_field_name_in_series_matrix, condition);
vector_X_axis{2} = extract_mean_vector_of_condition(stdized_table_X, H1_seriesTable,...
                                                condition_field_name_in_series_matrix, condition);
vector_Y_axis{1} = extract_mean_vector_of_condition(origin_tableY, M1_seriesTable,...
                                                condition_field_name_in_series_matrix, condition);
vector_Y_axis{2} = extract_mean_vector_of_condition(stdized_table_Y, M1_seriesTable,...
                                                condition_field_name_in_series_matrix, condition);
fig = figure;
%correlation between repeats
x=(0:1:15); y=(0:1:15);

subplot(1,2,1);

scatter(vector_X_axis{1},vector_Y_axis{1},'.');
xlim([0 15]); ylim([0 15]);
hold on
plot(x,y,'r');
set(gca, 'ytick',0:4:15, 'xtick',0:4:15, 'yticklabel', 0:4:15,...
    'xticklabel', 0:4:15, 'FontSize',10);

title('No Normalization')

subplot(1,2,2);

scatter(vector_X_axis{2},vector_Y_axis{2},'.');
xlim([0 15]); ylim([0 15]);
hold on
plot(x,y,'r');
set(gca, 'ytick',0:4:15, 'xtick',0:4:15, 'yticklabel', 0:4:15,...
    'xticklabel', 0:4:15, 'FontSize',10);

title('Mean Standartization')
saveas(fig,'mean_scatter_plots','epsc')          
% X_axis_mean_vecs = create_mean_vectors_for_condition_and_datasets(X_axis_exptable_name,...
%                             X_axis_series, condition_field_name_in_series_matrix, condition);
% Y_axis_mean_vecs = create_mean_vectors_for_condition_and_datasets(Y_axis_exptable_names,...
%                             Y_axis_series, condition_field_name_in_series_matrix, condition);




%                         
% datasets_titles = [{'H1-M1'}, {'H2-M2'}, {'H1-M2'}, {'H2-M1'}, {'H1-M2'}, {'M1-M2'}];
% normalization_titles = [{'No Normalization'}, {'Baseline'}, {'EB'}, {'DWD'}, {'XPN'}];
% 
% 
% [rows, columns] = size(X_axis_mean_vecs);                  
                        
% fig = figure
%     %correlation between repeats
% x=(0:1:15); y=(0:1:15);
% place=1; %each subplot position
% for i=1:rows
%     for j=1:columns
% 
%             subplot(rows,columns,place);
% 
%             scatter(X_axis_mean_vecs{i,j},Y_axis_mean_vecs{i,j},'.');
%             xlim([0 15]); ylim([0 15]);
%             hold on
%             plot(x,y,'r');
%             set(gca, 'ytick',0:4:15, 'xtick',0:4:15, 'yticklabel', 0:4:15,...
%                 'xticklabel', 0:4:15, 'FontSize',10);
%             
%             if i == 1
%             title(datasets_titles(j))    
%             end
%             place=place+1;
%             
%     end
% end

% saveas(fig,'scatter_plots','epsc')

function mean_vectors_cell = create_mean_vectors_for_condition_and_datasets(exptable_pre_norm_names,...
                            seriesTables, condition_field_name_in_series_matrix, condition)
    table_names = [exptable_pre_norm_names; strrep(exptable_pre_norm_names,'pre_norm','after_baseline');...
        strrep(exptable_pre_norm_names,'pre_norm','after_EB'); strrep(exptable_pre_norm_names,'pre_norm','after_DWD');...
        strrep(exptable_pre_norm_names,'pre_norm','afterXPN')];
    [rows_num, columns_num] = size(table_names);
    mean_vectors_cell = cell(rows_num, columns_num);
    for i = 1:rows_num
        for j = 1:columns_num
            table_name = table_names(i,j);
            table =  readtable(char(table_name{:}));
            mean_vectors_cell{i, j} = extract_mean_vector_of_condition(table, seriesTables{j},...
                                                    condition_field_name_in_series_matrix, condition);
        end
    end
end



function mean_vector_for_condition = extract_mean_vector_of_condition(exptable, seriesTable,...
                                                condition_field_name_in_series_matrix, condition)
    samples_name = filterByValues(seriesTable,seriesTable(:,1),...
                                  condition_field_name_in_series_matrix,condition);
    condition_samples = table2array(exptable(:,samples_name(1,:)));
    mean_vector_for_condition = mean(condition_samples,2);
end

function cellArray = readcell(path)
    T = readtable(path);
    cellArray = table2cell(T);
end