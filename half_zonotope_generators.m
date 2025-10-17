function reduced = half_zonotope_generators(fr, cutoff)
% 结合阈值筛选并缩减可行域生成元
G = fr.generators;
if isempty(G)
    reduced = G;
    return
end
importance = sum(abs(G), 1);
mask = importance >= cutoff.value;
if ~any(mask)
    [~, idx] = max(importance);
    mask(idx) = true;
end
reduced = G(:, mask);
end
