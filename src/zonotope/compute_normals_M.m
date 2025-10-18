function M = compute_normals_M(A_or_edges, b)
%COMPUTE_NORMALS_M Derive sorted unique unit normals for polygon/half-space.
%   M = COMPUTE_NORMALS_M(A, b) uses the rows of the half-space matrix A as
%   candidate outward normals (ignoring b except for size consistency) and
%   returns a 2-by-n matrix whose columns are unique, unit-norm, sorted by
%   polar angle.
%
%   M = COMPUTE_NORMALS_M(edges) accepts either polygon vertices (N-by-2)
%   or explicit edge vectors (2-by-K). In the vertex case, edges are
%   constructed by forward differencing the closed polygon. All candidate
%   normals are rotated by +90 degrees to point outward assuming clockwise
%   ordering.
%
%   Example
%   -------
%   A = [1 0; -1 0; 0 1; 0 -1]; b = [1; 0; 1; 0];
%   M = compute_normals_M(A, b);
%
%   See also REGION_TO_HALFSPACE

arguments
    A_or_edges (:,:) double
    b (:,1) double = []
end

if ~isempty(b)
    assert(size(A_or_edges,1) == numel(b), ...
        'Rows of A must match length of b.');
    normals = A_or_edges;
else
    data = A_or_edges;
    if size(data,2) == 2
        if any(abs(data(1,:) - data(end,:)) > 1e-9)
            data = [data; data(1,:)];
        end
        edges = diff(data,1,1);
        normals = zeros(size(edges));
        for k = 1:size(edges,1)
            edge = edges(k,:).';
            if norm(edge) <= eps
                continue;
            end
            n = [-edge(2); edge(1)];
            normals(k,:) = (n / norm(n)).';
        end
    elseif size(data,1) == 2
        edges = data;
        normals = zeros(size(edges));
        for k = 1:size(edges,2)
            edge = edges(:,k);
            if norm(edge) <= eps
                continue;
            end
            n = [-edge(2); edge(1)];
            normals(:,k) = n / norm(n);
        end
        normals = normals.';
    else
        error('Unsupported input dimensions for normals computation.');
    end
end

normals = normals(any(abs(normals) > eps, 2), :);
if isempty(normals)
    M = zeros(2,0);
    return;
end

% Normalise rows
norms = vecnorm(normals, 2, 2);
valid = norms > eps;
normals = normals(valid, :) ./ norms(valid);

angles = mod(atan2(normals(:,2), normals(:,1)), 2*pi);
[angles_sorted, order] = sort(angles);
normals = normals(order, :);

% Deduplicate by angle within tolerance
angle_tol = 1e-6;
unique_idx = [true; diff(angles_sorted) > angle_tol];
normals = normals(unique_idx, :);

M = normals.';
end
