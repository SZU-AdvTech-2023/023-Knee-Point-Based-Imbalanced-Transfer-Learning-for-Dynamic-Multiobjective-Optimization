index = [1,2, 3, 4, 5,6, 7, 9, 13]; % 指定要处理的索引数组

for i = 1:numel(index)
    subplot(4, 3, i); % 使用当前循环的索引更新 subplot
    
    num = index(i); % 获取当前迭代的索引值
    
    data = load(['MIGD-DF', num2str(num), '-P.txt']);
    
    % 绘制折线图
    hold on;
    plot(data, 'LineWidth', 1.5); % 绘制折线图，并设置线宽为 1.5
    
    % 使用 scatter 函数显示数据点，并设置颜色
    scatter(1:length(data), data, 50, 'b', 'filled'); % 'r' 代表红色
    
    title(sprintf('MIGD-DF%d-P', num)); % 设置标题
    xlabel('P'); % 设置 x 轴标签
    ylabel('MIGD'); % 设置 y 轴标签
    grid on; % 添加网格线
    %legend('Location', 'best'); % 添加图例
    hold off;

end
