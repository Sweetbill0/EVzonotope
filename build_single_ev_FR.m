function fr = build_single_ev_FR(sample, cfg)
% 基于单条样本生成近似的可行域多胞体
vector = sample(:)';
center = mean(vector);
deviations = vector - center;
segment_length = max(1, floor(numel(vector) / cfg.max_generators));
raw_generators = zeros(segment_length, cfg.max_generators);
used = 0;
for k = 1:cfg.max_generators
    idx_start = (k - 1) * segment_length + 1;
    idx_end = min(idx_start + segment_length - 1, numel(deviations));
    if idx_start > numel(deviations)
        break
    end
    block = deviations(idx_start:idx_end);
    padded = zeros(segment_length, 1);
    padded(1:numel(block)) = block(:);
    raw_generators(:, k) = padded;
    used = used + 1;
    if idx_end == numel(deviations)
        break
    end
end
generators = raw_generators(:, 1:used);
radius = norm(deviations, cfg.similarity_norm);
fr = struct();
fr.center = center;
fr.generators = generators;
fr.radius = radius;
fr.generator_length = segment_length;
end
