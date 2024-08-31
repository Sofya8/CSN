function build_folders_for_outputs()

    if ~exist('estimation_output','dir')
        mkdir('estimation_output')
    end

    if ~exist('estimation_output\filtered_lists','dir')
        mkdir('estimation_output\filtered_lists')
    end

    if ~exist('estimation_output\intersections','dir')
        mkdir('estimation_output\intersections')
    end

    if ~exist('estimation_output\separated_intersections','dir')
        mkdir('estimation_output\separated_intersections')
    end
end