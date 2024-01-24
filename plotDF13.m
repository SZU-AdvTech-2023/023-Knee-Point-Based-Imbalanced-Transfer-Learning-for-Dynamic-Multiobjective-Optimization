% % 第一个子图
% subplot(1, 2, 1); % 1 行 2 列，当前是第 1 个子图
% hold on;
% for i = 20
%     data = load(['KT-DF13-nt5-taut10environment', num2str(i), '-POF.txt']);
%     scatter3(data(:, 1), data(:, 2), data(:, 3), 'DisplayName', ['文件 ', num2str(i)]); % 使用前三列作为 x、y、z 轴
% end
% 
% % 添加图例
% legend('Location', 'best'); % 根据需要更改图例位置
% 
% % 添加标题和标签
% title('第一个子图'); % 更改图表标题
% xlabel('X 轴'); % 更改 x 轴标签
% ylabel('Y 轴'); % 更改 y 轴标签
% zlabel('Z 轴'); % 更改 z 轴标签
% 
% view(3); % 设置视角为三维
% hold off;

% 第二个子图
for num=7
    subplot(4, 4, num); 
    hold on;
    for i = 1:20
        data = load(['./Benchmark/pof/POF-nt10-taut5-DF',num2str(num) ,'-', num2str(i), '.txt']);
        %data = load(['KT-DF13-nt5-taut10environment', num2str(num), '-POF.txt']);
        if size(data,2) == 3
            scatter3(data(:, 1), data(:, 2), data(:, 3), 'DisplayName', ['文件 ', num2str(num)]); % 使用前三列作为 x、y、z 轴
        else
            plot(data(:,1), data(:,2), 'g');
            %scatter(data(:, 1), data(:, 2), 'DisplayName', ['文件 ', num2str(num)]); % 使用前三列作为 x、y、z 轴
        end
    end
    
    % 添加图例
    %legend('Location', 'best'); % 根据需要更改图例位置
    
    % 添加标题和标签
    title(sprintf('DF-%d', num)); % 更改图表标题
    xlabel('X 轴'); % 更改 x 轴标签
    ylabel('Y 轴'); % 更改 y 轴标签
    
    % 如果是三维数据，添加 z 轴标签
    if size(data, 2) == 3
        zlabel('Z 轴'); % 更改 z 轴标签
        view(3); % 设置视角为三维
    end

    hold off;
end
