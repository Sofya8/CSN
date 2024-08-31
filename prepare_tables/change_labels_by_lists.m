function target_table = change_labels_by_lists(target_table, labels_table, change_from_idx, to_idx )
    old_labels = target_table(:,1);
    lia = ismember (labels_table(:,change_from_idx), old_labels);
    labels_table = labels_table(lia,:);
    lia = ismember (old_labels, labels_table(:,change_from_idx));
    old_labels = old_labels(lia,:);
    target_table = target_table(lia,:);
    [~, locb] = ismember (labels_table(:,change_from_idx), old_labels);
    target_table = target_table(locb,:);
    target_table(:,1) = labels_table(:,to_idx);
end