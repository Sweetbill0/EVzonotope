function tests = tests()
%TESTS Unit tests for zonotope approximation pipeline.
%   Run with RUNTESTS to ensure core functionality remains correct.

tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(genpath('src'));
testCase.TestData.cfg = config();
end

function testRegionToHalfspaceContainsVertices(testCase)
poly = [0 0; 2 0; 2 1; 0 1; 0 0];
[A, b] = region_to_halfspace(poly);
vals = A * poly(1:end-1, :).';
verifyLessThanOrEqual(testCase, vals, b + 1e-9);

pt = mean(poly(1:4,:), 1).';
verifyLessThanOrEqual(testCase, A * pt, b + 1e-9);
end

function testWidthAlongNormalsMatchesProjection(testCase)
poly = [0 0; 1 0; 1 1; 0 1; 0 0];
[A, b] = region_to_halfspace(poly);
M = compute_normals_M(A, b);
DF = width_along_normals(A, b, M, testCase.TestData.cfg.linprog_opts);

for k = 1:size(M,2)
    dir = M(:,k);
    proj = poly(1:end-1,:) * dir;
    expected = max(proj) - min(proj);
    verifyEqual(testCase, DF(k), expected, 'AbsTol', 1e-6);
end
end

function testFitZonotopeLPLzConsistency(testCase)
cfg = testCase.TestData.cfg;
ev = struct('Cap_min', 5, 'Cap_max', 20, 'Cap_in', 10, 'Cap_exp', 18, ...
    'Pu', 6, 'Pd', 2, 't_in', 3, 't_out', 8);
poly = build_single_ev_region(ev, cfg.tgrid);
[A, b] = region_to_halfspace(poly);
M = compute_normals_M(A, b);
DF = width_along_normals(A, b, M, cfg.linprog_opts);
valid = ~isnan(DF);
M = M(:, valid);
DF = DF(valid);
Gu = build_generators_Gu(ev.Pu, ev.Pd, true);
[c, beta, report] = fit_zonotope_lp(A, b, Gu, M, DF, cfg.linprog_opts, cfg.norm); %#ok<ASGLU>
expected_Lz = 2 * abs(M.' * Gu) * beta;
verifyEqual(testCase, report.Lz, expected_Lz, 'AbsTol', 1e-6);
expected_Delta = mean(expected_Lz ./ DF);
verifyEqual(testCase, report.DeltaZ, expected_Delta, 'AbsTol', 1e-6);
end

function testMinkowskiSumDimensions(testCase)
z1 = struct('c', [1; 2], 'G', [1 0; 0 1], 'beta', [0.5; 0.4]);
z2 = struct('c', [0.5; 1], 'G', [0.2 0; 0 0.3], 'beta', [1; 1]);
agg = minkowski_sum_zonotopes([z1, z2]);
verifyEqual(testCase, agg.c, z1.c + z2.c, 'AbsTol', 1e-12);
verifySize(testCase, agg.G, [2, size(z1.G,2) + size(z2.G,2)]);
verifyEqual(testCase, numel(agg.beta), numel(z1.beta) + numel(z2.beta));
end

function testSamplingContainment(testCase)
A = [1 0; -1 0; 0 1; 0 -1];
b = [1; 0; 1; 0];
Gu = [0.4 0; 0 0.4];
c = [0.5; 0.5];
beta = [1; 1];
pts = sample_points_in_zonotope(c, Gu, beta, 1000);
inside = is_in_polytope(A, b, pts, 1e-6);
verifyGreaterThanOrEqual(testCase, mean(inside), 0.95);
verifyLessThanOrEqual(testCase, mean(inside), 1.0);
end
