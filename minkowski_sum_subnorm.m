function summary = minkowski_sum_subnorm(generators_a, generators_b, cfg)
% 计算生成元的米氏和并评估子范数差异
if isempty(generators_a)
    generators_a = zeros(size(generators_b, 1), 0);
end
if isempty(generators_b)
    generators_b = zeros(size(generators_a, 1), 0);
end
sum_generators = [generators_a, generators_b];
subnorm = norm(sum(generators_a, 2) - sum(generators_b, 2), cfg.similarity_norm);
summary = struct();
summary.generators = sum_generators;
summary.subnorm = subnorm;
end
