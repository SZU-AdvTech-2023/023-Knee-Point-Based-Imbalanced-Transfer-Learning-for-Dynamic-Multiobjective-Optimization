hold on;

for num = 1:13
    data = load(['./result/MIGD-DF', num2str(num), '-Svm.txt']);
    plot(data, 'LineWidth', 1.5)
    h1 = scatter(num, data, 50, 'r', 'filled'); % 'r' 代表红色
end

for num = 1:13
    data = load(['./result/MIGD-DF', num2str(num), '-SvmAndTree.txt']);
    plot(data, 'LineWidth', 1.5)
    h2 = scatter(num, data, 50, 'b', 'filled'); % 'r' 代表红色
end


legend([h1(1),h2(1)],'Svm', 'SvmAndTree','location', 'northwest');



% 添加标题和标签
title('Comparison of SVM and SVM&TREE'); % 设置图表标题
xlabel('DF'); % 设置 x 轴标签
ylabel('MIGD'); % 设置 y 轴标签
hold off;
