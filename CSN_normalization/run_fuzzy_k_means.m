function index = run_fuzzy_k_means(data, clasters_number)
    [~,U] = fcm(data,clasters_number);
    maxU = max(U);
    index = cell(clasters_number,1);
    for i=1:clasters_number
        index{i} = find(U(i,:) == maxU);
    end
end