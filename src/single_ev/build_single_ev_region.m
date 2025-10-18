function poly = build_single_ev_region(ev, tgrid)
%BUILD_SINGLE_EV_REGION Construct single-EV feasible region polygon.
%   poly = BUILD_SINGLE_EV_REGION(ev, tgrid) returns a clockwise-ordered set
%   of (time, energy) vertices describing the feasible "time-energy" region
%   bounded by the upper/lower envelopes defined in the specification.
%
%   Input
%   -----
%   ev : struct with fields Cap_min, Cap_max, Cap_in, Cap_exp, Pu, Pd, t_in,
%       t_out. All values use consistent units (hours for time, kWh for
%       energy, kW for power).
%   tgrid : row vector of window breakpoints used elsewhere in the
%       simulation (only its range [0, 24] is asserted here).
%
%   Output
%   ------
%   poly : N-by-2 array of vertices ordered clockwise (first column time,
%       second column energy). The polygon is suitable for conversion to
%       half-space form and downstream processing.
%
%   Example
%   -------
%   cfg = config();
%   ev = struct('Cap_min', 10, 'Cap_max', 40, 'Cap_in', 20, 'Cap_exp', 32, ...
%       'Pu', 7, 'Pd', 3, 't_in', 5, 't_out', 14);
%   poly = build_single_ev_region(ev, cfg.tgrid);
%
%   See also REGION_TO_HALFSPACE

arguments
    ev struct
    tgrid (1,:) double {mustBeNondecreasing}
end

must_have_fields = ["Cap_min","Cap_max","Cap_in","Cap_exp","Pu","Pd","t_in","t_out"];
assert(all(isfield(ev, must_have_fields)), 'Missing required EV fields.');

tol = 1e-6;
assert(ev.t_out - ev.t_in > -tol, 'Departure must not precede arrival.');
assert(ev.t_in >= -tol && ev.t_out <= 24 + tol, 'Times must stay within [0, 24].');
assert(ev.Cap_max - ev.Cap_min > tol, 'Battery span must be positive.');
assert(ev.Pu >= 0 && ev.Pd >= 0, 'Power limits must be non-negative.');
assert(ev.Cap_in >= 0 && ev.Cap_exp >= 0, 'Energies must be non-negative.');
assert(tgrid(1) <= ev.t_in + tol && tgrid(end) >= ev.t_out - tol, ...
    'tgrid must cover the activity window.');

% Short aliases
Cap_min = ev.Cap_min;
Cap_max = ev.Cap_max;
Cap_in = ev.Cap_in;
Cap_exp = min(max(ev.Cap_exp, Cap_min), Cap_max);
Pu = ev.Pu;
Pd = ev.Pd;
t_in = ev.t_in;
t_out = ev.t_out;

% Upper boundary construction ------------------------------------------------
upper_pts = [t_in, Cap_in];
if Pu <= tol
    % No ability to charge: energy stays at Cap_in, capped by Cap_max.
    upper_pts = unique_rows([upper_pts; t_out, min(Cap_in, Cap_max)]);
else
    t_c_max_end = t_in + max(0, (Cap_max - Cap_in)) / max(Pu, tol);
    if t_c_max_end > t_out + tol
        % Cannot reach Cap_max before departure.
        E_out = Cap_in + Pu * (t_out - t_in);
        upper_pts = unique_rows([upper_pts; t_out, min(E_out, Cap_max)]);
    else
        upper_pts = unique_rows([upper_pts; t_c_max_end, Cap_max; t_out, Cap_max]);
    end
end

% Lower boundary construction ------------------------------------------------
lower_pts = [t_in, Cap_in];
t_st = t_out - (Cap_exp - Cap_min) / max(Pu, tol);
t_force = clamp(t_st, t_in, t_out);

if Cap_in >= Cap_min - tol
    % Scenario I: arrival energy above minimum.
    if Pd <= tol
        % No discharge capability: hold then charge if needed.
        if t_force < t_out - tol
            lower_pts = unique_rows([lower_pts; t_force, Cap_in; t_out, Cap_exp]);
        else
            lower_pts = unique_rows([lower_pts; t_out, max(min(Cap_in, Cap_max), Cap_min)]);
        end
    else
        t_d_end = t_in + (Cap_in - Cap_min) / max(Pd, tol);
        Cap_out_d = Cap_in - Pd * (t_out - t_in);

        if t_force >= t_out - tol
            % Branch 1: forced charging starts at/after departure.
            if Cap_out_d <= Cap_min + tol
                lower_pts = unique_rows([lower_pts; t_d_end, Cap_min; t_out, Cap_min]);
            else
                E_out = max(Cap_min, Cap_out_d);
                lower_pts = unique_rows([lower_pts; t_out, E_out]);
            end
        elseif t_d_end <= t_force + tol
            % Branch 2: discharge hits Cap_min before forced charging begins.
            lower_pts = unique_rows([lower_pts; t_d_end, Cap_min]);
            if t_force > t_d_end + tol
                lower_pts = unique_rows([lower_pts; t_force, Cap_min]);
            end
            lower_pts = unique_rows([lower_pts; t_out, Cap_exp]);
        else
            % Branch 3: forced charging occurs before reaching Cap_min.
            E_force = Cap_in - Pd * (t_force - t_in);
            lower_pts = unique_rows([lower_pts; t_force, max(E_force, Cap_min); t_out, Cap_exp]);
        end
    end
else
    % Scenario II: arrival energy below minimum, immediate charging required.
    if Pu <= tol
        lower_pts = unique_rows([lower_pts; t_out, Cap_in]);
    else
        t_c_min_end = t_in + (Cap_min - Cap_in) / max(Pu, tol);
        if t_c_min_end >= t_out - tol
            E_out = Cap_in + Pu * (t_out - t_in);
            lower_pts = unique_rows([lower_pts; t_out, max(min(E_out, Cap_max), Cap_exp)]);
        elseif t_c_min_end <= t_force + tol
            % Line A: reach Cap_min before forced charging.
            lower_pts = unique_rows([lower_pts; t_c_min_end, Cap_min]);
            if t_force > t_c_min_end + tol
                lower_pts = unique_rows([lower_pts; t_force, Cap_min]);
            end
            lower_pts = unique_rows([lower_pts; t_out, Cap_exp]);
        else
            % Line B: forced charging starts before Cap_min is reached.
            lower_pts = unique_rows([lower_pts; t_c_min_end, Cap_min; t_out, Cap_exp]);
        end
    end
end

% Combine upper and lower boundaries (ensure clockwise order)
upper_pts = remove_consecutive_duplicates(upper_pts, tol);
lower_pts = remove_consecutive_duplicates(lower_pts, tol);

assert(abs(upper_pts(1,1) - lower_pts(1,1)) <= tol, 'Start times mismatch.');
assert(abs(upper_pts(end,1) - lower_pts(end,1)) <= tol, 'End times mismatch.');

poly = [upper_pts; flipud(lower_pts(2:end-1))];
poly = remove_consecutive_duplicates(poly, tol);

if any(abs(poly(1,:) - poly(end,:)) > tol)
    poly = [poly; poly(1,:)];
end

assert(size(poly,2) == 2, 'Polygon must be Nx2.');
assert(size(poly,1) >= 4, 'Polygon requires at least 4 vertices.');

end

function val = clamp(x, lo, hi)
val = min(max(x, lo), hi);
end

function arr = unique_rows(arr)
if isempty(arr)
    return;
end
[~, ia] = unique(round(arr, 9), 'rows', 'stable');
arr = arr(sort(ia), :);
end

function arr = remove_consecutive_duplicates(arr, tol)
if isempty(arr)
    return;
end
mask = [true; any(abs(diff(arr,1,1)) > tol, 2)];
arr = arr(mask, :);
end
