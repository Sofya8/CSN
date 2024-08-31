table_data = readtable('csc_index.xlsx');
b = figure;
%reordercats(categorical(table_data{1:6, 1}), table_data{1:6, 2:6})
bar(categorical(table_data{1:6, 1}), table_data{1:6, 2:6});
legend(table_data.Properties.VariableNames(2:6));
ylabel('CSC index');
saveas(b,'csc_index','epsc')

table_data = readtable('csc_index_cmaes.xlsx');
b = figure;
%reordercats(categorical(table_data{1:6, 1}), table_data{1:6, 2:6})
bar(categorical(table_data{1:6, 1}), table_data{1:6, 2:6});
legend(table_data.Properties.VariableNames(2:6));
ylabel('CSC index');
saveas(b,'csc_index_cmaes','epsc')