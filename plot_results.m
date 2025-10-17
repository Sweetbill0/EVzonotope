function fig = plot_results(report)
% 绘制管线的若干统计结果
entries = report.entries;
indices = cellfun(@(x) x.index, entries);
centers = cellfun(@(x) x.fr.center, entries);
deltas = cellfun(@(x) x.delta.score, entries);
subnorms = cellfun(@(x) x.aggregate_snapshot.subnorm, entries);
fig = figure('Name', 'EV Zonotope Summary', 'Visible', 'off');
subplot(3, 1, 1);
plot(indices, centers, '-o');
ylabel('中心');
subplot(3, 1, 2);
plot(indices, deltas, '-s');
ylabel('差值');
subplot(3, 1, 3);
plot(indices, subnorms, '-^');
ylabel('子范数');
xlabel('样本索引');
end
