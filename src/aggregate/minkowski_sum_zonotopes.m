function agg = minkowski_sum_zonotopes(zonotopes)
%MINKOWSKI_SUM_ZONOTOPES Aggregate a list of zonotopes via Minkowski sum.
%   agg = MINKOWSKI_SUM_ZONOTOPES(zonotopes) returns the Minkowski sum of
%   all input zonotopes. The input may be a struct array or cell array; each
%   element must contain fields 'c', 'G', and 'beta'. The output struct uses
%   the same field names.
%
%   Input
%   -----
%   zonotopes : Collection (struct array or cell array) of zonotope
%               descriptors. Each descriptor must include:
%                 c    - 2-by-1 centre vector
%                 G    - 2-by-d generator matrix
%                 beta - d-by-1 non-negative scales
%
%   Output
%   ------
%   agg : Struct with fields 'c', 'G', and 'beta' representing the summed
%         zonotope: centres add directly, generators concatenate, and beta
%         stacks vertically.
%
%   Example
%   -------
%   z1 = struct('c', [1; 2], 'G', eye(2), 'beta', [0.5; 0.3]);
%   z2 = struct('c', [3; 1], 'G', eye(2), 'beta', [0.2; 0.4]);
%   agg = minkowski_sum_zonotopes([z1, z2]);
%
%   See also FIT_ZONOTOPE_LP

arguments
    zonotopes
end

if isempty(zonotopes)
    agg = struct('c', zeros(2,1), 'G', zeros(2,0), 'beta', zeros(0,1));
    return;
end

if iscell(zonotopes)
    zonos = [zonotopes{:}];
elseif isstruct(zonotopes)
    zonos = zonotopes;
else
    error('zonotopes must be a struct array or cell array.');
end

required_fields = {'c','G','beta'};
for k = 1:numel(zonos)
    assert(all(isfield(zonos(k), required_fields)), ...
        'Each zonotope must contain fields c, G, beta.');
    validateattributes(zonos(k).c, {'double'}, {'size', [2,1], 'finite'});
    validateattributes(zonos(k).G, {'double'}, {'finite'});
    validateattributes(zonos(k).beta, {'double'}, {'column', 'nonnegative'});
    if isempty(zonos(k).G)
        assert(isempty(zonos(k).beta), 'Empty generator matrix requires empty beta.');
    else
        assert(size(zonos(k).G,1) == 2, 'Generator matrix must have two rows.');
        assert(size(zonos(k).G,2) == numel(zonos(k).beta), ...
            'Mismatch between generator columns and beta length.');
    end
end

c_sum = zeros(2,1);
G_blocks = {};
beta_blocks = {};

for k = 1:numel(zonos)
    c_sum = c_sum + zonos(k).c;
    if ~isempty(zonos(k).G)
        G_blocks{end+1} = zonos(k).G; %#ok<AGROW>
        beta_blocks{end+1} = zonos(k).beta; %#ok<AGROW>
    end
end

if isempty(G_blocks)
    G_sum = zeros(2,0);
    beta_sum = zeros(0,1);
else
    G_sum = [G_blocks{:}];
    beta_sum = vertcat(beta_blocks{:});
end

agg = struct('c', c_sum, 'G', G_sum, 'beta', beta_sum);
end
