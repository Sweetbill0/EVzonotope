function tests()
% 简易测试覆盖核心函数
cfg = config();
[dataset, ~] = simulate_ev_dataset(cfg);
sample = dataset.fetch(1);
fr = build_single_ev_FR(sample, cfg);
assert(isfield(fr, 'generators'));
halfspace = fit_to_halfspace(fr, cfg);
assert(isfield(halfspace, 'normal'));
cutoff = compute_cutoff_of_normals(halfspace, cfg);
reduced = half_zonotope_generators(fr, cutoff);
summary = minkowski_sum_subnorm(fr.generators, reduced, cfg);
assert(isfield(summary, 'generators'));
delta = compute_similarity_delta_2(fr, fr, cfg);
assert(delta.score >= 0);
report = aggregator_pipeline();
assert(numel(report.entries) == cfg.sample_size);
disp('tests completed');
end
