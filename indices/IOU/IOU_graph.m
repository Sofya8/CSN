table_data = readtable('IOU_index.xlsx');
f = figure;
b = bar(categorical(table_data{1:6, 1}), table_data{1:6, 2:9});
for i=1:4
    b(2*i-1).LineWidth = 0.0001;
    b(2*i).LineWidth = 2;
    b(2*i).FaceColor = [1 1 1];
    b(2*i).EdgeColor = b(2*i-1).FaceColor ;
end
legend(table_data.Properties.VariableNames(2:9));
ylabel('IOU index');
saveas(f,'IOU_rate','epsc')

table_data = readtable('IOU_index_cmaes.xlsx');
f = figure;
b = bar(categorical(table_data{1:6, 1}), table_data{1:6, 2:9});
for i=1:4
    b(2*i-1).LineWidth = 0.0001;
    b(2*i).LineWidth = 2;
    b(2*i).FaceColor = [1 1 1];
    b(2*i).EdgeColor = b(2*i-1).FaceColor ;
end
legend(table_data.Properties.VariableNames(2:9));
ylabel('IOU index');
saveas(f,'IOU_rate_cmaes','epsc')

table_data = readtable('average_IOU_index_cmaes.xlsx');
b = figure;
bar(categorical(table_data{1:6, 1}), table_data{1:6, 2:5});
legend(table_data.Properties.VariableNames(2:5));
ylabel('Average IOU index');
saveas(b,'average_IOU_rate','epsc')
