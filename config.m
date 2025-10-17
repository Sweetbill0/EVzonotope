function cfg = config()
% 配置结构体涵盖仿真范围与算法参数
cfg = struct();
cfg.random_seed = 2025;
cfg.dataset_rows = 30629;
cfg.dataset_cols = 20160;
cfg.significant_digits = 7;
cfg.max_generators = 64;
cfg.reduction_tolerance = 1e-3;
cfg.similarity_norm = 2;
cfg.sample_size = 32;
cfg.linprog_options = optimoptions('linprog', 'Display', 'none');
end
