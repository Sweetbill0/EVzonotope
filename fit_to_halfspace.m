function halfspace = fit_to_halfspace(fr, cfg)
% 通过线性规划为可行域寻找紧致半空间
G = fr.generators;
if isempty(G)
    normal = ones(1, 1);
    offset = fr.center;
    exitflag = 1;
else
    dim = size(G, 1);
    cols = size(G, 2);
    f = [zeros(dim, 1); 1];
    A = [-G', ones(cols, 1); G', ones(cols, 1)];
    b = fr.radius * ones(2 * cols, 1);
    lb = [-inf(dim, 1); 0];
    [solution, ~, exitflag] = linprog(f, A, b, [], [], lb, [], [], cfg.linprog_options);
    if exitflag <= 0 || isempty(solution)
        normal = zeros(dim, 1);
        normal(1) = 1;
        offset = fr.center + fr.radius;
    else
        normal = solution(1:dim);
        offset = fr.center + solution(end);
    end
end
halfspace = struct();
halfspace.normal = normal;
halfspace.offset = offset;
halfspace.exitflag = exitflag;
end
