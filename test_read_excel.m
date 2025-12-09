% 测试Excel读取函数（只读取S11和S21，不含频率）
% 使用方法：
%   1. 将您的Excel文件路径填入下面
%   2. 设置频率范围参数
%   3. 运行此脚本
%   4. 查看读取结果和绘图

%% 配置
excel_file = 'your_file.xlsx';  % 修改为您的Excel文件名

% 可选参数
options.sheet = 1;           % 工作表编号
options.start_row = 1;       % 数据起始行（无标题时从1开始）

% 频率设置（需要手动提供）
freq_start = 8e9;    % 起始频率 (Hz)，例如8 GHz
freq_stop = 12e9;    % 终止频率 (Hz)，例如12 GHz
% 数据点数会根据Excel行数自动确定

%% 读取数据
fprintf('========== 测试Excel读取 ==========\n\n');

[S11, S21] = read_excel_s_params(excel_file, options);

% 生成频率向量（等间隔）
n_points = length(S11);
freq = linspace(freq_start, freq_stop, n_points)';

fprintf('\n生成频率向量: %.2f ~ %.2f GHz (%d点)\n', ...
        freq_start/1e9, freq_stop/1e9, n_points);

%% 显示前5个数据点
fprintf('\n前5个数据点:\n');
fprintf('%-12s %-25s %-25s\n', '频率(GHz)', 'S11', 'S21');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:min(5, length(freq))
    fprintf('%-12.4f %-25s %-25s\n', ...
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
fprintf('接下来可以使用诊断函数进行材料参数计算。\n\n');

%% 提示：如何使用诊断函数
fprintf('========== 下一步操作 ==========\n');
fprintf('如需计算材料参数，先创建临时S2P文件：\n\n');
fprintf('  %% 保存为S2P格式\n');
fprintf('  write_s2p_file(''temp.s2p'', freq, S11, S21);\n\n');
fprintf('  %% 然后运行诊断函数\n');
fprintf('  results = calculate_s2p_diagnostic(''temp.s2p'', 1.4e-3);\n\n');
