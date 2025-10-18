function [users, records] = simulate_ev_dataset(cfg)
%SIMULATE_EV_DATASET Generate synthetic EV population and windowed records.
% 用途: 基于配置生成车辆参数并展开到所有时间窗。
% 输入: cfg 结构体，由 CONFIG 返回。
% 输出: users（结构数组）与 records（表），分别为用户级参数与时间窗记录。
% 示例: cfg = config(); [users, records] = simulate_ev_dataset(cfg);
%   head(records)
%
%   See also CONFIG

arguments
    cfg (1,1) struct
end

validateattributes(cfg.N_users, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.tgrid, {'numeric'}, {'row','nondecreasing'});
assert(cfg.tgrid(1) == 0 && cfg.tgrid(end) == 24, ...
    'Time grid must cover [0, 24].');
assert(cfg.Nwin == numel(cfg.tgrid) - 1, ...
    'cfg.Nwin inconsistent with cfg.tgrid.');

rng(cfg.seed);

segNames = cfg.seg_names;
numSeg = numel(segNames);
weights = zeros(1, numSeg);
for k = 1:numSeg
    weights(k) = cfg.seg.(segNames{k}).weight;
end
weights = weights ./ sum(weights);
segmentIndex = randsample(numSeg, cfg.N_users, true, weights);

users = repmat(struct( ...
    'id', 0, ...
    'segment', "", ...
    'Cap_max', 0, ...
    'Cap_min', 0, ...
    'Cap_in', 0, ...
    'Cap_exp', 0, ...
    'Pu', 0, ...
    'Pd', 0, ...
    't_in', 0, ...
    'T_in', 0, ...
    't_out', 0 ...
    ), cfg.N_users, 1);

for i = 1:cfg.N_users
    segName = segNames{segmentIndex(i)};
    params = cfg.seg.(segName);

    Cap_max = uniform_range(params.Cap_max_rng);
    Cap_min = params.Cap_min_ratio * Cap_max;
    Cap_in = Cap_min + (Cap_max - Cap_min) * betarnd(params.Cap_in_beta(1), params.Cap_in_beta(2));
    ratio = params.Cap_exp_rule.ratio;
    tol = params.Cap_exp_rule.tolerance;
    Cap_exp = ratio * Cap_max * (1 + tol * (2*rand() - 1));
    Cap_exp = min(max(Cap_exp, Cap_min), Cap_max);

    Pu = uniform_range(params.Pu_rng);
    if rand() <= params.V2G_rate
        Pd = uniform_range(params.Pd_rng);
    else
        Pd = 0;
    end

    t_in = sample_truncated_mog(params.t_in_mog_params, 0, 24);
    T_in = sample_duration(params.T_in_rng);
    t_out = min(t_in + T_in, 24);
    T_in = t_out - t_in;

    users(i).id = i;
    users(i).segment = segName;
    users(i).Cap_max = Cap_max;
    users(i).Cap_min = Cap_min;
    users(i).Cap_in = Cap_in;
    users(i).Cap_exp = Cap_exp;
    users(i).Pu = Pu;
    users(i).Pd = Pd;
    users(i).t_in = t_in;
    users(i).T_in = T_in;
    users(i).t_out = t_out;
end

windowStarts = cfg.tgrid(1:end-1);
windowEnds = cfg.tgrid(2:end);
Nwin = cfg.Nwin;

userIdx = repelem((1:cfg.N_users)', Nwin);
winIdx = repmat((1:Nwin)', cfg.N_users, 1);
winStart = repmat(windowStarts, cfg.N_users, 1);
winStart = winStart(:);
winEnd = repmat(windowEnds, cfg.N_users, 1);
winEnd = winEnd(:);

Cap_max_all = kron([users.Cap_max]', ones(Nwin,1));
Cap_min_all = kron([users.Cap_min]', ones(Nwin,1));
Cap_in_all = kron([users.Cap_in]', ones(Nwin,1));
Cap_exp_all = kron([users.Cap_exp]', ones(Nwin,1));
Pu_all = kron([users.Pu]', ones(Nwin,1));
Pd_all = kron([users.Pd]', ones(Nwin,1));
t_in_all = kron([users.t_in]', ones(Nwin,1));
t_out_all = kron([users.t_out]', ones(Nwin,1));
segment_all = strings(cfg.N_users * Nwin, 1);
for i = 1:cfg.N_users
    idxRange = (i-1)*Nwin + (1:Nwin);
    segment_all(idxRange) = users(i).segment;
end

is_active = (winEnd > t_in_all) & (winStart < t_out_all);

records = table(userIdx, segment_all, winIdx, winStart, winEnd, ...
    t_in_all, t_out_all, is_active, Cap_max_all, Cap_min_all, ...
    Cap_in_all, Cap_exp_all, Pu_all, Pd_all, ...
    'VariableNames', {'user_id','segment','window_id','t_start','t_end', ...
    't_in','t_out','is_active','Cap_max','Cap_min','Cap_in','Cap_exp','Pu','Pd'});

end

function x = uniform_range(bounds)
validateattributes(bounds, {'numeric'}, {'vector','numel',2});
lo = bounds(1);
hi = bounds(2);
assert(hi > lo, 'Invalid range bounds.');
x = lo + (hi - lo) * rand();
end

function t = sample_truncated_mog(params, lo, hi)
weights = [params.weight];
weights = weights(:)';
weights = weights / sum(weights);
mu = [params.mu];
sigma = [params.sigma];
valid = false;
while ~valid
    component = randsample(numel(weights), 1, true, weights);
    sample = mu(component) + sigma(component) * randn();
    if sample >= lo && sample <= hi
        valid = true;
    end
end
t = sample;
end

function d = sample_duration(desc)
switch desc.type
    case 'uniform'
        d = uniform_range(desc.bounds);
    case 'mixture_uniform'
        weights = desc.weights(:)';
        weights = weights / sum(weights);
        component = randsample(numel(weights), 1, true, weights);
        bounds = desc.bounds(component, :);
        d = uniform_range(bounds);
    otherwise
        error('Unsupported duration descriptor type: %s', desc.type);
end
end
