function results = run_demo_single()
%RUN_DEMO_SINGLE Reproduce single-EV zonotope fitting workflow demonstration.
%   results = RUN_DEMO_SINGLE() executes the single-vehicle pipeline: load
%   configuration, sample a small set of EVs, construct feasible regions,
%   fit baseline and improved zonotopes, compute quality metrics and
%   generate comparison plots saved to ./out/figs.
%
%   Output
%   ------
%   results : struct containing per-vehicle diagnostics for the baseline
%       and improved generator sets, including DeltaZ similarity scores,
%       containment rates, and area ratios.
%
%   Example
%   -------
%   out = run_demo_single();
%   disp(out.summary);
%
%   See also CONFIG, SIMULATE_EV_DATASET, BUILD_SINGLE_EV_REGION,
%   REGION_TO_HALFSPACE, BUILD_GENERATORS_GU, FIT_ZONOTOPE_LP

arguments
end

addpath(genpath('src'));

cfg = config();
rng(cfg.seed, 'twister');

[users, ~] = simulate_ev_dataset(cfg);
num_candidates = numel(users);
assert(num_candidates >= 3, 'Dataset must contain at least three users.');

num_select = randi([1, 3]);
perm = randperm(num_candidates, num_select);
selected_users = users(perm);

baseline_metrics = repmat(initialize_metric_struct(), num_select, 1);
improved_metrics = repmat(initialize_metric_struct(), num_select, 1);

if ~exist(fullfile('out','figs'), 'dir')
    mkdir(fullfile('out','figs'));
end

for i = 1:num_select
    ev = selected_users(i);
    poly = build_single_ev_region(ev, cfg.tgrid);
    [A, b] = region_to_halfspace(poly);

    M = compute_normals_M(A, b);
    DF = width_along_normals(A, b, M, cfg.linprog_opts);
    valid_idx = ~isnan(DF);
    if ~any(valid_idx)
        error('No valid directional widths obtained for user %d.', ev.id);
    end
    M = M(:, valid_idx);
    DF = DF(valid_idx);

    [baseline_metrics(i), baseline_zono] = fit_and_evaluate(ev, poly, A, b, M, DF, false, cfg);
    [improved_metrics(i), improved_zono] = fit_and_evaluate(ev, poly, A, b, M, DF, true, cfg);

    fig = figure('Color', 'w', 'Visible', 'off');
    ax = axes('Parent', fig);
    hold(ax, 'on');
    plot_region_polytope(ax, poly);
    plot_zonotope(ax, baseline_zono.c, baseline_zono.G, baseline_zono.beta, struct( ...
        'LineColor', [0.15, 0.70, 0.25], 'LineWidth', 1.5, 'LineStyle', '--'));
    plot_zonotope(ax, improved_zono.c, improved_zono.G, improved_zono.beta, struct( ...
        'LineColor', [0.85, 0.33, 0.10], 'LineWidth', 2.0, 'LineStyle', '-'));
    legend(ax, {'Feasible region', 'Baseline zonotope', 'Improved zonotope'}, 'Location', 'best');
    title(ax, sprintf('User %d (%s)', ev.id, ev.segment));
    hold(ax, 'off');

    save_path = fullfile('out', 'figs', sprintf('single_demo_user_%04d.png', ev.id));
    exportgraphics(fig, save_path, 'Resolution', 150);
    close(fig);
end

summary = struct();
summary.num_users = num_select;
summary.baseline = summarize_metrics(baseline_metrics);
summary.improved = summarize_metrics(improved_metrics);

print_summary(summary);

results = struct('baseline', {baseline_metrics}, 'improved', {improved_metrics}, 'summary', summary);

end

function metric = initialize_metric_struct()
metric = struct('DeltaZ', NaN, 'containment', NaN, 'area_ratio', NaN, ...
    'center', [], 'G', [], 'beta', [], 'report', []);
end

function [metric, zono] = fit_and_evaluate(ev, poly, A, b, M, DF, use_improved, cfg)
Gu = build_generators_Gu(ev.Pu, ev.Pd, use_improved);

[c, beta, report] = fit_zonotope_lp(A, b, Gu, M, DF, cfg.linprog_opts, cfg.norm);

num_samples = 1000;
samples = sample_points_in_zonotope(c, Gu, beta, num_samples);
inside = is_in_polytope(A, b, samples, cfg.norm.eps);
containment = mean(inside);

poly_area = polygon_area(poly);
zono_area = zonotope_area(c, Gu, beta);
area_ratio = zono_area / poly_area;

metric = struct('DeltaZ', report.DeltaZ, 'containment', containment, ...
    'area_ratio', area_ratio, 'center', c, 'G', Gu, 'beta', beta, 'report', report);

zono = struct('c', c, 'G', Gu, 'beta', beta);

end

function area = polygon_area(poly)
if norm(poly(1,:) - poly(end,:), 2) > 1e-9
    poly = [poly; poly(1,:)];
end
area = polyarea(poly(:,1), poly(:,2));
end

function area = zonotope_area(c, G, beta) %#ok<INUSD>
d = size(G,2);
if d == 0
    area = 0;
    return;
end
scaled_G = G .* beta.';
num_vertices = 2^d;
vertices = zeros(num_vertices, 2);
for idx = 0:num_vertices-1
    signs = 2 * bitget(idx, 1:d).' - 1;
    point = c + scaled_G * signs;
    vertices(idx+1, :) = point.';
end
[~, area] = convhull(vertices(:,1), vertices(:,2)); %#ok<ASGLU>
end

function summary = summarize_metrics(metrics)
DeltaZ = [metrics.DeltaZ]';
containment = [metrics.containment]';
area_ratio = [metrics.area_ratio]';

summary = struct();
summary.DeltaZ = describe_vector(DeltaZ);
summary.containment = describe_vector(containment);
summary.area_ratio = describe_vector(area_ratio);
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

function print_summary(summary)
fprintf('Single-EV demo across %d users\n', summary.num_users);
print_stats('Baseline   ΔZ', summary.baseline.DeltaZ);
print_stats('Improved   ΔZ', summary.improved.DeltaZ);
print_stats('Baseline   containment', summary.baseline.containment);
print_stats('Improved   containment', summary.improved.containment);
print_stats('Baseline   area ratio', summary.baseline.area_ratio);
print_stats('Improved   area ratio', summary.improved.area_ratio);
end

function print_stats(label, stats)
fprintf('  %-22s mean=%6.3f  [P25=%6.3f  P50=%6.3f  P75=%6.3f]\n', ...
    label, stats.mean, stats.p25, stats.p50, stats.p75);
end
