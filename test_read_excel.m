% 测试Excel读取函数
% 使用方法：
%   1. 将您的Excel文件路径填入下面
%   2. 运行此脚本
%   3. 查看读取结果和绘图

%% 配置
excel_file = 'your_file.xlsx';  % 修改为您的Excel文件名

% 可选参数
options.sheet = 1;           % 工作表编号
options.start_row = 2;       % 数据起始行（跳过标题）
options.freq_unit = 'GHz';   % 频率单位

%% 读取数据
fprintf('========== 测试Excel读取 ==========\n\n');

[freq, S11, S21] = read_excel_s_params(excel_file, options);

%% 显示前5个数据点
fprintf('\n前5个数据点:\n');
fprintf('%-12s %-20s %-20s\n', '频率(GHz)', 'S11', 'S21');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:min(5, length(freq))
    fprintf('%-12.4f %-20s %-20s\n', ...
            freq(i)/1e9, ...
            sprintf('%.4f%+.4fi', real(S11(i)), imag(S11(i))), ...
            sprintf('%.4f%+.4fi', real(S21(i)), imag(S21(i))));
end

%% 绘制S参数
figure('Name', 'Excel数据读取验证', 'Position', [100, 100, 1200, 600]);

% S11幅度和相位
subplot(2, 2, 1);
plot(freq/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (GHz)');
ylabel('幅度 (dB)');
title('S11幅度');

subplot(2, 2, 3);
plot(freq/1e9, angle(S11)*180/pi, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (GHz)');
ylabel('相位 (度)');
title('S11相位');

% S21幅度和相位
subplot(2, 2, 2);
plot(freq/1e9, 20*log10(abs(S21)), 'r-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (GHz)');
ylabel('幅度 (dB)');
title('S21幅度');

subplot(2, 2, 4);
plot(freq/1e9, angle(S21)*180/pi, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('频率 (GHz)');
ylabel('相位 (度)');
title('S21相位');

fprintf('\n如果图形显示正常，说明读取成功！\n');
fprintf('接下来可以使用这些数据进行材料参数计算。\n\n');
