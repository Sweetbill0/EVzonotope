function h = plot_region_polytope(ax, poly, opts)
%PLOT_REGION_POLYTOPE Visualise feasible region polygon.
%   h = PLOT_REGION_POLYTOPE(ax, poly) renders the polygon specified by
%   vertices POLY onto the axes AX using a semi-transparent patch.
%
%   h = PLOT_REGION_POLYTOPE(ax, poly, opts) allows additional styling via
%   the structure OPTS with fields:
%       FaceColor  : RGB triplet (default [0.20 0.45 0.80])
%       FaceAlpha  : Transparency in [0, 1] (default 0.25)
%       EdgeColor  : RGB triplet (default [0.10 0.25 0.55])
%       LineWidth  : Positive scalar edge width (default 1.5)
%
%   Input
%   -----
%   ax   : Target axes handle. The function asserts AX is valid and ready
%          for plotting.
%   poly : N-by-2 array of (time, energy) vertices ordered clockwise or
%          counter-clockwise. The polygon need not be closed; closure is
%          enforced internally.
%   opts : (optional) struct controlling patch appearance.
%
%   Output
%   ------
%   h    : Patch graphics object handle representing the polygon.
%
%   Example
%   -------
%   cfg = config();
%   ev = struct('Cap_min', 10, 'Cap_max', 40, 'Cap_in', 20, 'Cap_exp', 30, ...
%       'Pu', 7, 'Pd', 3, 't_in', 4, 't_out', 12);
%   poly = build_single_ev_region(ev, cfg.tgrid);
%   fig = figure('Visible', 'off');
%   ax = axes('Parent', fig);
%   h = plot_region_polytope(ax, poly);
%
%   See also BUILD_SINGLE_EV_REGION, PLOT_ZONOTOPE

arguments
    ax (1,1) matlab.graphics.axis.Axes
    poly (:,2) double
    opts.FaceColor (1,3) double {mustBeNonnegative, mustBeLessThanOrEqual(opts.FaceColor, 1)} = [0.20, 0.45, 0.80]
    opts.FaceAlpha (1,1) double {mustBeGreaterThanOrEqual(opts.FaceAlpha, 0), mustBeLessThanOrEqual(opts.FaceAlpha, 1)} = 0.25
    opts.EdgeColor (1,3) double {mustBeNonnegative, mustBeLessThanOrEqual(opts.EdgeColor, 1)} = [0.10, 0.25, 0.55]
    opts.LineWidth (1,1) double {mustBeNonnegative} = 1.5
end

validateattributes(poly, {'double'}, {'ncols', 2, 'finite', 'real'});
assert(size(poly,1) >= 3, 'Polygon must contain at least three vertices.');

if any(isnan(poly), 'all') || any(isinf(poly), 'all')
    error('Polygon coordinates must be finite.');
end

tol = 1e-9;
poly = ensure_closed_polygon(poly, tol);

hold_state = ishold(ax);
hold(ax, 'on');

h = patch('Parent', ax, ...
    'XData', poly(:,1), ...
    'YData', poly(:,2), ...
    'FaceColor', opts.FaceColor, ...
    'FaceAlpha', opts.FaceAlpha, ...
    'EdgeColor', opts.EdgeColor, ...
    'LineWidth', opts.LineWidth);

grid(ax, 'on');
axis(ax, 'tight');
ax.Layer = 'top';
ax.Box = 'on';
ax.XLabel.String = 'Time [h]';
ax.YLabel.String = 'Energy [kWh]';

if ~hold_state
    hold(ax, 'off');
end

end

function poly = ensure_closed_polygon(poly, tol)
if norm(poly(1,:) - poly(end,:), 2) > tol
    poly = [poly; poly(1,:)];
end

mask = [true; any(abs(diff(poly,1,1)) > tol, 2)];
poly = poly(mask, :);

if size(poly,1) < 3
    error('Polygon degenerates after duplicate removal.');
end
end
