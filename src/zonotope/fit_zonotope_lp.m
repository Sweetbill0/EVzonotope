function [c, beta, report] = fit_zonotope_lp(A, b, Gu, M, DF, opts, normSpec)
%FIT_ZONOTOPE_LP Fit inner zonotope via linear programming.
%   [c, beta, report] = FIT_ZONOTOPE_LP(A, b, Gu, M, DF, opts) solves a
%   linear program to maximise the average directional similarity between
%   the polytope defined by Ax <= b and the zonotope Z(c, {beta_j g_j}). The
%   decision variables comprise the centre c (2-by-1) and the non-negative
%   scaling vector beta.
%
%   The problem is solved in a normalised [0,1]^2 coordinate system for
%   numerical stability and mapped back to the original units on return.
%
%   Input
%   -----
%   A, b : Half-space representation of the polygon (A*x <= b).
%   Gu   : 2-by-d generator matrix.
%   M    : 2-by-nf matrix of unit normal directions.
%   DF   : nf-by-1 vector of polytope widths along columns of M.
%   opts : optimoptions structure for LINPROG (optional).
%   normSpec : struct with fields t_min, E_min, scale_t, scale_E (optional).
%
%   Output
%   ------
%   c     : 2-by-1 zonotope centre in original units.
%   beta  : d-by-1 non-negative generator scales.
%   report: struct containing LP diagnostics and similarity metrics.
%
%   Example
%   -------
%   cfg = config();
%   ev = struct('Cap_min', 10, 'Cap_max', 40, 'Cap_in', 20, 'Cap_exp', 32, ...
%       'Pu', 7, 'Pd', 3, 't_in', 5, 't_out', 14);
%   poly = build_single_ev_region(ev, cfg.tgrid);
%   [A, b] = region_to_halfspace(poly);
%   Gu = build_generators_Gu(ev.Pu, ev.Pd, true);
%   M = compute_normals_M(A, b);
%   DF = width_along_normals(A, b, M, cfg.linprog_opts);
%   [c, beta, report] = fit_zonotope_lp(A, b, Gu, M, DF, cfg.linprog_opts, cfg.norm);
%
%   See also BUILD_GENERATORS_GU, WIDTH_ALONG_NORMALS, LINPROG

arguments
    A (:,:) double
    b (:,1) double
    Gu (2,:) double
    M (2,:) double
    DF (:,1) double
    opts = optimoptions('linprog', 'Algorithm', 'dual-simplex', 'Display', 'none')
    normSpec struct = struct('t_min', 0, 'E_min', 0, 'scale_t', 1, 'scale_E', 1)
end

assert(size(A,2) == 2, 'This implementation expects 2-D polytopes.');
assert(size(M,1) == 2, 'Normal matrix must have two rows.');
assert(numel(DF) == size(M,2), 'DF length must match number of normals.');

d = size(Gu,2);
assert(d >= 1, 'Generator matrix must contain at least one column.');

nf = size(M,2);
valid_idx = ~isnan(DF) & DF > 0;
assert(any(valid_idx), 'At least one valid direction is required.');
M = M(:, valid_idx);
DF = DF(valid_idx);
nf = numel(DF);

% Normalisation -------------------------------------------------------------
D = diag([normSpec.scale_t; normSpec.scale_E]);
D_inv = diag(1 ./ diag(D));
offset = [normSpec.t_min; normSpec.E_min];

A_n = A * D;
b_n = b - A * offset;
G_n = D_inv * Gu;
M_n = D.' * M;

% Renormalise normals to unit length and adjust DF accordingly
norms = vecnorm(M_n, 2, 1);
valid_norms = norms > eps;
if ~any(valid_norms)
    error('fit_zonotope_lp:degenerateNormals', 'No valid normals after normalisation.');
end
M_n = M_n(:, valid_norms);
M = M(:, valid_norms);
DF_orig = DF(valid_norms);
DF = DF_orig ./ norms(valid_norms).';
nf = numel(DF);
M_n = M_n ./ norms(valid_norms);

% Objective coefficients ----------------------------------------------------
abs_AG = abs(A_n * G_n);
abs_MG = abs(M_n.' * G_n);
weights = (2 / nf) * (1 ./ DF);
beta_coeff = abs_MG.' * weights;

f = zeros(2 + d, 1);
f(3:end) = -beta_coeff;

% Inequality constraints ----------------------------------------------------
A_ineq = [A_n, abs_AG];
b_ineq = b_n;

lb = [-inf(2,1); zeros(d,1)];

[x_opt, fval, exitflag, output] = linprog(f, A_ineq, b_ineq, [], [], lb, [], [], opts);

if exitflag <= 0
    warning('fit_zonotope_lp:linprogFailed', 'linprog did not converge: %s', output.message);
end

c_scaled = x_opt(1:2);
beta = x_opt(3:end);

c = D * c_scaled + offset;

Lz = 2 * abs(M.' * Gu) * beta;
DeltaZ = mean(Lz ./ DF_orig);

report = struct();
report.exitflag = exitflag;
report.fval = fval;
report.Lz = Lz;
report.DeltaZ = DeltaZ;
report.DF = DF_orig;
report.linprog_output = output;
report.normalised = struct('c', c_scaled, 'A', A_n, 'b', b_n, 'G', G_n);

end
