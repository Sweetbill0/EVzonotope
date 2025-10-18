function fig = plot_agg_demo(agg_zono, segment_stats, save_path, agg_poly)
%PLOT_AGG_DEMO Visualise aggregated zonotope and segment contributions.
%   fig = PLOT_AGG_DEMO(agg_zono, segment_stats) creates a two-panel figure
%   comparing the aggregated zonotope footprint against an optional
%   reference polygon and displaying per-segment contribution statistics.
%
%   fig = PLOT_AGG_DEMO(agg_zono, segment_stats, save_path) additionally
%   writes the figure to SAVE_PATH (PNG). When SAVE_PATH is empty, the
%   figure remains unsaved.
%
%   fig = PLOT_AGG_DEMO(agg_zono, segment_stats, save_path, agg_poly)
%   overlays the polygon AGG_POLY (if provided) in the first subplot using
%   PLOT_REGION_POLYTOPE.
%
%   Input
%   -----
%   agg_zono      : struct with fields 'c', 'G', 'beta' representing the
%                   aggregated zonotope parameters.
%   segment_stats : table containing a categorical or string 'segment'
%                   column and a numeric contribution column. Recognised
%                   column names for contributions are 'contribution',
%                   'value', or 'count'.
%   save_path     : (optional) string or char path for saving the figure.
%                   Defaults to 'out/figs/aggregate_demo.png'.
%   agg_poly      : (optional) N-by-2 array of polygon vertices for the
%                   aggregated feasible region.
%
%   Output
%   ------
%   fig           : Figure handle containing the generated plots.
%
%   Example
%   -------
%   agg_zono = struct('c', [10; 200], 'G', [1 0; 0 1], 'beta', [2; 1]);
%   stats = table(["Private"; "Taxi"], [120; 80], ...
%       'VariableNames', {'segment','contribution'});
%   fig = plot_agg_demo(agg_zono, stats);
%
%   See also PLOT_REGION_POLYTOPE, PLOT_ZONOTOPE

arguments
    agg_zono struct
    segment_stats table
    save_path = fullfile('out', 'figs', 'aggregate_demo.png')
    agg_poly (:,2) double = []
end

required_fields = {'c', 'G', 'beta'};
assert(all(isfield(agg_zono, required_fields)), 'agg_zono requires fields c, G, beta.');
validateattributes(agg_zono.c, {'double'}, {'size', [2, 1], 'finite', 'real'});
validateattributes(agg_zono.G, {'double'}, {'2d', 'finite', 'real'});
assert(size(agg_zono.G,1) == 2, 'agg_zono.G must have exactly two rows.');
validateattributes(agg_zono.beta, {'double'}, {'column', 'nonnegative', 'finite'});
assert(numel(agg_zono.beta) == size(agg_zono.G,2), ...
    'Length of beta must match number of generator columns.');

assert(istable(segment_stats), 'segment_stats must be a table.');
assert(any(strcmpi(segment_stats.Properties.VariableNames, 'segment')), ...
    'segment_stats must include a ''segment'' column.');
segments = segment_stats.segment;

if iscell(segments)
    segments = string(segments);
elseif iscategorical(segments)
    segments = string(segments);
elseif ~isstring(segments)
    error('segment column must be convertible to string.');
end

value_col = get_value_column(segment_stats);
values = segment_stats.(value_col);
validateattributes(values, {'double'}, {'vector', 'real', 'finite'});
assert(numel(values) == numel(segments), 'segment/value size mismatch.');
values = values(:);
segments = segments(:);

fig = figure('Color', 'w', 'Visible', 'off');
fig.Position(3:4) = [900, 420];

ax1 = subplot(1,2,1, 'Parent', fig);
legend_handles = gobjects(0);
legend_labels = strings(0,1);
if ~isempty(agg_poly)
    h_poly = plot_region_polytope(ax1, agg_poly, struct('FaceColor', [0.30, 0.60, 0.90], ...
        'FaceAlpha', 0.20, 'EdgeColor', [0.20, 0.40, 0.70], 'LineWidth', 1.2));
    legend_handles(end+1) = h_poly; %#ok<AGROW>
    legend_labels(end+1) = "Aggregated region"; %#ok<AGROW>
end

h_zono = plot_zonotope(ax1, agg_zono.c, agg_zono.G, agg_zono.beta, struct( ...
    'LineColor', [0.85, 0.33, 0.10], 'LineWidth', 2.2, 'LineStyle', '-', ...
    'NumAngles', 512));
legend_handles(end+1) = h_zono; %#ok<AGROW>
legend_labels(end+1) = "Aggregated zonotope"; %#ok<AGROW>
ax1.Title.String = sprintf('Aggregated Zonotope (|G|=%d)', size(agg_zono.G,2));
legend(ax1, legend_handles, legend_labels, 'Location', 'best');

ax2 = subplot(1,2,2, 'Parent', fig);
[segments_sorted, order] = sort(segments);
values_sorted = values(order);
bar(ax2, categorical(segments_sorted), values_sorted, 'FaceColor', [0.20, 0.60, 0.50]);
ax2.Title.String = 'Segment Contributions';
value_label = char(value_col);
if ~isempty(value_label)
    value_label = lower(value_label);
    value_label(1) = upper(value_label(1));
else
    value_label = 'Contribution';
end
ax2.YLabel.String = value_label;
ax2.XLabel.String = 'Segment';
ax2.YGrid = 'on';

if ~exist('out/figs', 'dir')
    mkdir('out/figs');
end

if ~(isstring(save_path) || ischar(save_path))
    error('save_path must be a string or character vector.');
end

if strlength(string(save_path)) > 0
    exportgraphics(fig, save_path, 'Resolution', 150);
end

end

function col = get_value_column(tbl)
candidates = {'contribution', 'value', 'count'};
for k = 1:numel(candidates)
    matches = strcmpi(tbl.Properties.VariableNames, candidates{k});
    if any(matches)
        col = tbl.Properties.VariableNames{find(matches, 1, 'first')};
        return;
    end
end
error('segment_stats must contain one of the columns: contribution, value, count.');
end
