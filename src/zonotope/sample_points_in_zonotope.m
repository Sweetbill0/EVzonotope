function samples = sample_points_in_zonotope(c, G, beta, num_samples)
%SAMPLE_POINTS_IN_ZONOTOPE Uniformly sample points from a zonotope.
%   samples = SAMPLE_POINTS_IN_ZONOTOPE(c, G, beta, num_samples) draws
%   NUM_SAMPLES random points from the zonotope Z(c, {beta_j g_j}) by
%   sampling independent coefficients uniformly in [-1, 1].
%
%   Input
%   -----
%   c          : 2-by-1 centre of the zonotope.
%   G          : 2-by-d generator matrix.
%   beta       : d-by-1 non-negative scaling vector.
%   num_samples: scalar positive integer specifying the number of points.
%
%   Output
%   ------
%   samples : 2-by-num_samples matrix; each column is a sampled point.
%
%   Example
%   -------
%   c = [5; 30]; G = [1 0; 0 1]; beta = [0.5; 0.4];
%   pts = sample_points_in_zonotope(c, G, beta, 1000);
%
%   See also IS_IN_POLYTOPE

arguments
    c (2,1) double {mustBeFinite, mustBeReal}
    G (2,:) double {mustBeFinite, mustBeReal}
    beta (:,1) double {mustBeNonnegative, mustBeFinite}
    num_samples (1,1) double {mustBeInteger, mustBePositive}
end

d = size(G,2);
assert(numel(beta) == d, 'Length of beta must match generator count.');

scaled_G = G .* beta.';
Xi = -1 + 2 * rand(d, num_samples);
samples = c + scaled_G * Xi;
end
