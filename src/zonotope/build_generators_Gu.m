function Gu = build_generators_Gu(Pu, Pd, use_improved_generators)
%BUILD_GENERATORS_GU Construct generator matrix for 2-D zonotopes.
%   Gu = BUILD_GENERATORS_GU(Pu, Pd, use_improved_generators) returns the
%   2-by-d generator matrix for the time-energy zonotope. The baseline
%   configuration uses the canonical basis vectors (time-only and
%   energy-only directions). The improved configuration augments these with
%   charging and discharging slope directions normalised to unit length.
%
%   Input
%   -----
%   Pu : scalar non-negative charging power limit (kW).
%   Pd : scalar non-negative discharging power limit (kW). Zero denotes
%        absence of V2G capability.
%   use_improved_generators : logical flag; true adds slope-aligned
%        directions, false returns the baseline {e1, e2}.
%
%   Output
%   ------
%   Gu : 2-by-d matrix whose columns are generator directions. Columns are
%        ordered as [time-axis, energy-axis, charge-slope, discharge-slope?].
%
%   Example
%   -------
%   G0 = build_generators_Gu(7, 0, false);   % baseline rectangle
%   G1 = build_generators_Gu(7, 3, true);    % augmented with slope vectors
%
%   See also FIT_ZONOTOPE_LP

arguments
    Pu (1,1) double {mustBeFinite, mustBeNonnegative}
    Pd (1,1) double {mustBeFinite, mustBeNonnegative}
    use_improved_generators (1,1) logical = true
end

tol = 1e-9;

e1 = [1; 0];
e2 = [0; 1];
Gu = [e1, e2];

if use_improved_generators
    vc = normalise_vector([1; Pu]);

    if Pd > tol
        vd = normalise_vector([-1; Pd]);
        Gu = [Gu, vc, vd];
    else
        Gu = [Gu, vc];
    end

    Gu = unique_columns(Gu, tol);
end

end

function v = normalise_vector(v)
nrm = norm(v);
if nrm <= eps
    error('Generator direction must be non-zero.');
end
v = v / nrm;
end

function M = unique_columns(M, tol)
if isempty(M)
    return;
end
cols = size(M, 2);
keep = true(1, cols);
for i = 2:cols
    if ~keep(i)
        continue;
    end
    for j = 1:i-1
        if keep(j) && norm(M(:,i) - M(:,j)) <= tol
            keep(i) = false;
            break;
        end
    end
end
M = M(:, keep);
end
