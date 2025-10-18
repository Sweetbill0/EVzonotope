function DF = width_along_normals(A, b, M, opts)
%WIDTH_ALONG_NORMALS Compute polytope width along specified directions.
%   DF = WIDTH_ALONG_NORMALS(A, b, M) returns the directional widths of the
%   polytope defined by {x | A*x <= b} along each column of M. Each column
%   of M must be a (near-)unit vector representing an outward normal
%   direction. The width is computed as max m' x - min m' x via linear
%   programming.
%
%   DF = WIDTH_ALONG_NORMALS(A, b, M, opts) allows passing custom options
%   created with OPTIMOPTIONS for LINPROG. When omitted, dual-simplex with
%   suppressed output is used.
%
%   Input
%   -----
%   A : m-by-n constraint matrix.
%   b : m-by-1 vector defining the right-hand side.
%   M : n-by-k matrix of direction vectors (ideally unit-norm).
%   opts : (optional) optimisation options for LINPROG.
%
%   Output
%   ------
%   DF : k-by-1 vector of directional widths. Entries are NaN for
%       directions that are numerically degenerate (width < 1e-6) or if the
%       LP fails to converge.
%
%   Example
%   -------
%   A = [1 0; -1 0; 0 1; 0 -1];
%   b = [1; 0; 1; 0];
%   M = [1 0; 0 1]';
%   DF = width_along_normals(A, b, M');
%
%   See also LINPROG

arguments
    A (:,:) double
    b (:,1) double
    M (:,:) double
    opts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'none')
end

[m, n] = size(A);
assert(numel(b) == m, 'Dimension mismatch between A and b.');
assert(size(M,1) == n, 'Direction matrix must have %d rows.', n);

k = size(M,2);
DF = nan(k,1);

fallback_opts = optimoptions('linprog', 'Algorithm', 'interior-point', 'Display', 'none');
tol = 1e-6;

for h = 1:k
    dir = M(:,h);
    nrm = norm(dir);
    assert(nrm > eps, 'Direction vectors must be non-zero.');
    if abs(nrm - 1) > tol
        dir = dir / nrm;
    end

    [~, fval_max, exitflag_max] = linprog(-dir, A, b, [], [], [], [], [], opts);
    if exitflag_max <= 0
        [~, fval_max, exitflag_max] = linprog(-dir, A, b, [], [], [], [], [], fallback_opts);
    end

    [~, fval_min, exitflag_min] = linprog(dir, A, b, [], [], [], [], [], opts);
    if exitflag_min <= 0
        [~, fval_min, exitflag_min] = linprog(dir, A, b, [], [], [], [], [], fallback_opts);
    end

    if exitflag_max <= 0 || exitflag_min <= 0 || isinf(fval_max) || isinf(fval_min)
        DF(h) = NaN;
        continue;
    end

    width = -fval_max - fval_min;
    if width < tol
        DF(h) = NaN;
    else
        DF(h) = width;
    end
end

end
