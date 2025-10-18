function inside = is_in_polytope(A, b, points, tol)
%IS_IN_POLYTOPE Check whether points satisfy polytope inequalities.
%   inside = IS_IN_POLYTOPE(A, b, points) returns a logical row vector
%   indicating whether each column of POINTS lies within the polytope
%   described by A*x <= b. POINTS may be 2-by-N or N-by-2; the function will
%   automatically orient the matrix appropriately.
%
%   inside = IS_IN_POLYTOPE(A, b, points, tol) allows specifying a
%   tolerance added to b when checking inequalities (default 1e-6).
%
%   Example
%   -------
%   A = [1 0; -1 0; 0 1; 0 -1];
%   b = [1; 0; 1; 0];
%   pts = [0.5 0.2 1.1; 0.5 0.1 0.8];
%   inside = is_in_polytope(A, b, pts);
%
%   See also SAMPLE_POINTS_IN_ZONOTOPE

arguments
    A (:,:) double
    b (:,1) double
    points (:,:) double
    tol (1,1) double {mustBeNonnegative} = 1e-6
end

if size(points,1) ~= size(A,2)
    if size(points,2) == size(A,2)
        points = points.';
    else
        error('Point dimensionality must match columns of A.');
    end
end

values = A * points;
inside = all(values <= b + tol, 1);
end
