function h = plot_zonotope(ax, center, G, beta, opts)
%PLOT_ZONOTOPE Draw a 2-D zonotope footprint.
%   h = PLOT_ZONOTOPE(ax, center, G, beta) renders the boundary of the
%   zonotope Z(center, {beta_k g_k}) on the axes AX as a line plot.
%
%   h = PLOT_ZONOTOPE(ax, center, G, beta, opts) enables name-value style
%   customisation via structure OPTS with fields:
%       LineColor : RGB triplet (default [0.85 0.33 0.10])
%       LineWidth : Positive scalar (default 2)
%       LineStyle : Line style string (default '-')
%       NumAngles : Integer >=16 controlling the angular resolution used to
%                   approximate the zonotope boundary (default 256).
%
%   Input
%   -----
%   ax     : Target axes handle.
%   center : 2-by-1 vector representing zonotope centre.
%   G      : 2-by-d generator matrix.
%   beta   : d-by-1 non-negative scaling vector.
%   opts   : Optional structure controlling plot appearance.
%
%   Output
%   ------
%   h      : Line object handle returned by PLOT.
%
%   Example
%   -------
%   ax = axes('Parent', figure('Visible', 'off'));
%   center = [5; 40];
%   G = [1 0 0.5; 0 1 0.5];
%   beta = [0.5; 0.4; 0.2];
%   h = plot_zonotope(ax, center, G, beta);
%
%   See also PLOT_REGION_POLYTOPE

arguments
    ax (1,1) matlab.graphics.axis.Axes
    center (2,1) double {mustBeFinite, mustBeReal}
    G (2,:) double {mustBeFinite, mustBeReal}
    beta (:,1) double {mustBeNonnegative, mustBeFinite}
    opts.LineColor (1,3) double {mustBeNonnegative, mustBeLessThanOrEqual(opts.LineColor, 1)} = [0.85, 0.33, 0.10]
    opts.LineWidth (1,1) double {mustBeNonnegative} = 2.0
    opts.LineStyle (1,:) char = '-'
    opts.NumAngles (1,1) double {mustBeInteger, mustBePositive} = 256
end

validateattributes(center, {'double'}, {'size', [2, 1]});
assert(size(G,1) == 2, 'Generator matrix must have two rows.');
d = size(G,2);
assert(numel(beta) == d, 'beta must match number of generator columns.');

assert(opts.NumAngles >= 16, 'NumAngles must be at least 16.');

if any(isnan(beta)) || any(isinf(beta))
    error('beta must contain finite values.');
end

vertices = compute_outline(center, G, beta, opts.NumAngles);

hold_state = ishold(ax);
hold(ax, 'on');

h = plot(ax, vertices(:,1), vertices(:,2), ...
    'Color', opts.LineColor, ...
    'LineWidth', opts.LineWidth, ...
    'LineStyle', opts.LineStyle);

if ~hold_state
    hold(ax, 'off');
end

end

function vertices = compute_outline(center, G, beta, num_angles)
d = size(G,2);
if d == 0
    vertices = [center.'; center.'];
    return;
end

angles = linspace(0, 2*pi, num_angles + 1);
angles(end) = [];
points = zeros(numel(angles), 2);

for k = 1:numel(angles)
    dir = [cos(angles(k)); sin(angles(k))];
    proj = G.' * dir;
    signs = sign(proj);
    signs(signs == 0) = 1;
    coeff = beta .* signs;
    point = center + G * coeff;
    points(k, :) = point.';
end

% Remove potential duplicates and compute convex hull for ordering.
[~, ia] = unique(round(points, 9), 'rows', 'stable');
points = points(sort(ia), :);

if size(points,1) < 3
    vertices = [points; points(1,:)];
    return;
end

[~, hull_idx] = convhull(points(:,1), points(:,2));
vertices = points(hull_idx, :);
end
