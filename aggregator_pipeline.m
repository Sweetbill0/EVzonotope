function report = aggregator_pipeline()
% 串联数据生成与多面体分析的流水线
cfg = config();
[dataset, meta] = simulate_ev_dataset(cfg);
indices = round(linspace(1, dataset.rows, cfg.sample_size));
entries = cell(numel(indices), 1);
aggregate = [];
for idx = 1:numel(indices)
    sample_index = indices(idx);
    sample = dataset.fetch(sample_index);
    fr = build_single_ev_FR(sample, cfg);
    halfspace = fit_to_halfspace(fr, cfg);
    cutoff = compute_cutoff_of_normals(halfspace, cfg);
    reduced = half_zonotope_generators(fr, cutoff);
    if isempty(aggregate)
        aggregate = struct();
        aggregate.fr = fr;
        aggregate.halfspace = halfspace;
        aggregate.cutoff = cutoff;
        aggregate.reduced = reduced;
        aggregate.subnorm = 0;
        aggregate.weight = 1;
    else
        combination = minkowski_sum_subnorm(aggregate.fr.generators, reduced, cfg);
        aggregate.fr.generators = combination.generators;
        aggregate.fr.center = (aggregate.fr.center * aggregate.weight + fr.center) / (aggregate.weight + 1);
        aggregate.subnorm = combination.subnorm;
        aggregate.halfspace = halfspace;
        aggregate.cutoff = cutoff;
        aggregate.reduced = reduced;
        aggregate.weight = aggregate.weight + 1;
    end
    delta = compute_similarity_delta_2(aggregate.fr, fr, cfg);
    entry = struct();
    entry.index = sample_index;
    entry.fr = fr;
    entry.halfspace = halfspace;
    entry.cutoff = cutoff;
    entry.reduced = reduced;
    entry.delta = delta;
    entry.aggregate_snapshot = aggregate;
    entries{idx} = entry;
end
report = struct();
report.cfg = cfg;
report.meta = meta;
report.dataset = dataset;
report.entries = entries;
report.aggregate = aggregate;
end
