function [A, b] = region_to_halfspace(poly)
%REGION_TO_HALFSPACE Convert polygon vertices to half-space representation.
%   [A, b] = REGION_TO_HALFSPACE(poly) computes outward unit normals for each
%   polygon edge and returns inequality constraints A*x <= b that describe
%   the convex polygon. The input polygon is expected to be ordered
%   clockwise; if not, it is reoriented internally.
%
%   Input
%   -----
%   poly : N-by-2 array of polygon vertices in (time, energy) space. The
%       polygon must be closed (first vertex equals last) or will be closed
%       automatically.
%
%   Output
%   ------
%   A : M-by-2 matrix where each row is an outward unit normal vector of an
%       edge.
%   b : M-by-1 vector with offsets such that the polygon satisfies A*x <= b
%       exactly (up to numerical tolerance).
%
%   Example
%   -------
%   square = [0 0; 1 0; 1 1; 0 1; 0 0];
%   [A, b] = region_to_halfspace(square);
%   assert(all(A * [0.5; 0.5] <= b + 1e-12));
%
%   See also BUILD_SINGLE_EV_REGION

arguments
    poly (:,2) double
end

tol = 1e-9;

assert(size(poly,1) >= 3, 'Polygon must provide at least three vertices.');

poly = remove_duplicate_vertices(poly, tol);
if any(abs(poly(1,:) - poly(end,:)) > tol)
    poly = [poly; poly(1,:)];
end

if polygon_signed_area(poly) > 0
    poly = flipud(poly);
end

num_edges = size(poly,1) - 1;
A = zeros(num_edges, 2);
b = zeros(num_edges, 1);

for k = 1:num_edges
    p1 = poly(k, :).';
    p2 = poly(k+1, :).';
    edge = p2 - p1;
    if norm(edge) <= tol
        error('Consecutive polygon vertices must not coincide.');
    end
    normal = [-edge(2); edge(1)];
    norm_val = norm(normal);
    assert(norm_val > tol, 'Degenerate edge encountered.');
    normal = normal / norm_val;
    A(k, :) = normal.';
    b(k) = normal.' * p1;
end

if any(A * poly(1:end-1,:).' - b > 1e-6)
    error('Polygon vertices violate constructed half-space constraints.');
end

end

function poly = remove_duplicate_vertices(poly, tol)
mask = [true; any(abs(diff(poly,1,1)) > tol, 2)];
poly = poly(mask, :);
if size(poly,1) >= 2 && norm(poly(1,:) - poly(end,:), 2) <= tol
    poly = poly(1:end-1, :);
end
end

function area = polygon_signed_area(poly)
x = poly(:,1);
y = poly(:,2);
area = 0.5 * sum(x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1));
end
