function [dataset, meta] = simulate_ev_dataset(cfg)
% 构造电动车样本生成器并产生元数据
rng(cfg.random_seed, 'twister');
trend = linspace(-0.25, 0.75, cfg.dataset_cols);
bias_vector = sin((1:cfg.dataset_rows)' * pi / cfg.dataset_rows);
row_phase = rand(cfg.dataset_rows, 1);
col_phase = rand(1, cfg.dataset_cols);
col_scale = rand(1, cfg.dataset_cols);
meta = struct();
meta.trend = trend;
meta.bias = bias_vector;
meta.seed = cfg.random_seed;
meta.row_phase = row_phase;
meta.col_phase = col_phase;
meta.col_scale = col_scale;
    function row = build_row(index)
        modulation = bias_vector(index);
        base_pattern = sin(row_phase(index) + col_phase);
        wave_pattern = cos(modulation + col_scale);
        row = base_pattern + modulation * trend + wave_pattern .* modulation;
        row = round(row, cfg.significant_digits, 'significant');
    end
dataset = struct();
dataset.rows = cfg.dataset_rows;
dataset.cols = cfg.dataset_cols;
dataset.trend = trend;
dataset.bias = bias_vector;
dataset.fetch = @build_row;
dataset.row_phase = row_phase;
dataset.col_phase = col_phase;
dataset.col_scale = col_scale;
end
