function main()
% 执行完整的多面体仿真流程
report = aggregator_pipeline();
plot_results(report);
entry_count = numel(report.entries);
centers = cellfun(@(x) x.fr.center, report.entries);
mean_center = mean(centers);
subnorms = cellfun(@(x) x.aggregate_snapshot.subnorm, report.entries);
summary = struct();
summary.sample_count = entry_count;
summary.mean_center = mean_center;
summary.final_subnorm = subnorms(end);
disp(summary);
end
