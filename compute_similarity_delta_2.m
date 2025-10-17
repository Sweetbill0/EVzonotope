function delta = compute_similarity_delta_2(reference, candidate, cfg)
% 衡量两个区域之间的相似度差值
center_diff = abs(reference.center - candidate.center);
ref_sum = sum(reference.generators, 2);
cand_sum = sum(candidate.generators, 2);
if isempty(ref_sum)
    ref_sum = zeros(size(cand_sum));
end
if isempty(cand_sum)
    cand_sum = zeros(size(ref_sum));
end
generator_diff = norm(ref_sum - cand_sum, cfg.similarity_norm);
score = center_diff + generator_diff;
delta = struct();
delta.center = center_diff;
delta.generator = generator_diff;
delta.score = score;
end
