function new_matrix = move_blocks(matrix, index, alphas, betas, clusters_number)
    
    [rows, cols] = size(matrix);
    new_matrix = zeros(rows, cols);
    for i = 1:clusters_number
        block_index = index{i};
        new_matrix(block_index,:) = alphas(i)*matrix(block_index,:)+betas(i);
    end
    
end