  for num=1:14
    subplot(4, 4, num);
    data = load(['KT-DF',num2str(num),'-nt5-taut10-IGD.txt']);
    
    % 绘制折线图
    plot(data, 'LineWidth', 1.5); % 绘制折线图，并设置线宽为 1.5
    hold on; % 保持图形，以便在其上添加 scatter 数据点
    
    % 使用 scatter 函数显示数据点，并设置颜色
    scatter(1:length(data), data, 50, 'b', 'filled'); % 'r' 代表红色
    
    title(sprintf('DF%d-nt5-taut10',num)); % 设置标题
    xlabel('n'); % 设置 x 轴标签
    ylabel('MIGD'); % 设置 y 轴标签
    grid on; % 添加网格线
    %legend('Location', 'best'); % 添加图例
    hold off; % 结束标记状态
  end

