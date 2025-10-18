function results = run_demo_aggregate()
%RUN_DEMO_AGGREGATE Aggregate zonotopes across time windows and segments.
%   results = RUN_DEMO_AGGREGATE() executes the full population workflow:
%   configure the simulation, generate all EV records, fit zonotopes for
%   each active user-window entry, perform segment-wise Minkowski sums per
%   time window, and collect similarity statistics.
%
%   Output
%   ------
%   results : struct containing per-record zonotope parameters, aggregated
%       zonotopes per window/segment, summary statistics, and runtime
%       measurements.
%
%   Example
%   -------
%   out = run_demo_aggregate();
%   disp(out.summary.DeltaZ);
%
%   See also CONFIG, SIMULATE_EV_DATASET, FIT_ZONOTOPE_LP, MINKOWSKI_SUM_ZONOTOPES

arguments
end

addpath(genpath('src'));

cfg = config();
rng(cfg.seed, 'twister');

[users, records] = simulate_ev_dataset(cfg); %#ok<NASGU>
active_records = records(records.is_active, :);
n_active = height(active_records);

fprintf('Total active records: %d\n', n_active);

zono_entries = repmat(struct('window_id', 0, 'segment', "", 'zono', [], ...
    'DeltaZ', NaN, 'report', []), n_active, 1);

fit_timer = tic;
for i = 1:n_active
    rec = active_records(i, :);
    ev = struct('Cap_min', rec.Cap_min, 'Cap_max', rec.Cap_max, ...
        'Cap_in', rec.Cap_in, 'Cap_exp', rec.Cap_exp, 'Pu', rec.Pu, ...
        'Pd', rec.Pd, 't_in', rec.t_in, 't_out', rec.t_out);

    poly = build_single_ev_region(ev, cfg.tgrid);
    [A, b] = region_to_halfspace(poly);
    M = compute_normals_M(A, b);
    DF = width_along_normals(A, b, M, cfg.linprog_opts);
    valid_idx = ~isnan(DF);
    if ~any(valid_idx)
        warning('Record %d produced no valid widths; skipping.', i);
        continue;
    end
    M = M(:, valid_idx);
    DF = DF(valid_idx);

    Gu = build_generators_Gu(rec.Pu, rec.Pd, cfg.use_improved_generators);
    [c, beta, report] = fit_zonotope_lp(A, b, Gu, M, DF, cfg.linprog_opts, cfg.norm);

    zono_entries(i).window_id = rec.window_id;
    zono_entries(i).segment = rec.segment;
    zono_entries(i).zono = struct('c', c, 'G', Gu, 'beta', beta);
    zono_entries(i).DeltaZ = report.DeltaZ;
    zono_entries(i).report = report;
end
fit_time = toc(fit_timer);

valid_mask = ~arrayfun(@(s) isempty(s.zono), zono_entries);
zono_entries = zono_entries(valid_mask);
active_records = active_records(valid_mask, :);
n_active = numel(zono_entries);

fprintf('Successful zonotopes: %d (fit time %.2f s)\n', n_active, fit_time);

window_ids = unique(active_records.window_id);
num_windows = numel(window_ids);

agg_per_window = repmat(struct('window_id', 0, 'segment_zonos', [], ...
    'total_zono', []), num_windows, 1);

for w = 1:num_windows
    wid = window_ids(w);
    idx_w = find(active_records.window_id == wid);
    segments = active_records.segment(idx_w);
    unique_segments = unique(segments);
    segment_zonos = repmat(struct('segment', "", 'count', 0, 'zono', []), numel(unique_segments), 1);
    segment_results = cell(numel(unique_segments), 1);

    for s = 1:numel(unique_segments)
        seg = unique_segments(s);
        idx_seg = idx_w(segments == seg);
        zonos = [zono_entries(idx_seg).zono];
        segment_zonos(s).segment = seg;
        segment_zonos(s).count = numel(idx_seg);
        if isempty(zonos)
            segment_zonos(s).zono = struct('c', zeros(2,1), 'G', zeros(2,0), 'beta', zeros(0,1));
        else
            segment_zonos(s).zono = minkowski_sum_zonotopes(zonos);
        end
        segment_results{s} = segment_zonos(s).zono;
    end

    total_zono = minkowski_sum_zonotopes(segment_results(~cellfun('isempty', segment_results)));
    agg_per_window(w).window_id = wid;
    agg_per_window(w).segment_zonos = segment_zonos;
    agg_per_window(w).total_zono = total_zono;
end

DeltaZ_vals = [zono_entries.DeltaZ]';
DeltaZ_stats = describe_vector(DeltaZ_vals);

summary = struct();
summary.num_records = n_active;
summary.DeltaZ = DeltaZ_stats;
summary.fit_time = fit_time;

results = struct();
results.entries = zono_entries;
results.agg_per_window = agg_per_window;
results.summary = summary;

% Prepare contribution table and figure for window with maximum records
[~, max_idx] = max(arrayfun(@(w) sum([w.segment_zonos.count]), agg_per_window));
if ~isempty(max_idx)
    chosen = agg_per_window(max_idx);
    contrib_segments = string({chosen.segment_zonos.segment});
    contrib_counts = [chosen.segment_zonos.count]';
    contrib_table = table(contrib_segments', contrib_counts, ...
        'VariableNames', {'segment','count'});
    fig = plot_agg_demo(chosen.total_zono, contrib_table, ...
        fullfile('out', 'figs', sprintf('aggregate_window_%02d.png', chosen.window_id)));
    close(fig);
end

fprintf('ΔZ mean=%.3f  [P25=%.3f, P50=%.3f, P75=%.3f]\n', ...
    summary.DeltaZ.mean, summary.DeltaZ.p25, summary.DeltaZ.p50, summary.DeltaZ.p75);

end

function stats = describe_vector(x)
x = x(:);
x = x(~isnan(x));
if isempty(x)
    stats = struct('mean', NaN, 'p25', NaN, 'p50', NaN, 'p75', NaN);
    return;
end
stats = struct('mean', mean(x), 'p25', prctile(x, 25), ...
    'p50', prctile(x, 50), 'p75', prctile(x, 75));
end
