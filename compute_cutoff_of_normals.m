function cutoff = compute_cutoff_of_normals(halfspace, cfg)
% 依据法向量幅值计算截断阈值
magnitudes = abs(halfspace.normal(:));
if isempty(magnitudes)
    cutoff_value = 0;
    limit = 0;
else
    sorted = sort(magnitudes, 'descend');
    limit = min(numel(sorted), cfg.max_generators);
    cutoff_value = sorted(limit);
end
cutoff = struct();
cutoff.value = cutoff_value;
cutoff.limit = limit;
cutoff.magnitudes = magnitudes;
end
