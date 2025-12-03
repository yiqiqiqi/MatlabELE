%% S2P文件材料参数计算 - 使用示例（X波段专用）
%
% 本脚本演示如何使用 calculate_from_s2p 函数
% 从S2P文件计算材料的相对磁导率、相对介电常数和损耗角正切
%
% 固定参数（X波段）:
%   - 截止波长 λc = 45.7 mm
%   - 波导宽边 a = 22.85 mm
%   - 适用频率: 8.2 - 12.4 GHz
%
% 唯一需要设置的参数: 样品厚度 d

clear all;
close all;
clc;

%% ========== 参数设置 ==========

% S2P文件路径（请修改为你的文件路径）
s2p_file = 'your_material_test.s2p';  % ← 修改这里

% 样品厚度 (单位: 米)
d = 2e-3;  % 例如: 2mm = 2e-3 m  ← 修改这里

%% ========== 执行计算 ==========

fprintf('╔════════════════════════════════════════════════════╗\n');
fprintf('║    X波段材料电磁参数计算（S2P文件）              ║\n');
fprintf('╚════════════════════════════════════════════════════╝\n\n');

fprintf('输入参数:\n');
fprintf('  S2P文件: %s\n', s2p_file);
fprintf('  样品厚度: %.2f mm\n\n', d*1000);

try
    % 调用计算函数
    results = calculate_from_s2p(s2p_file, d);

    %% ========== 显示部分结果 ==========
    fprintf('\n========== 详细结果（部分频率点）==========\n');

    % 显示5个均匀分布的频率点
    n_points = length(results.frequency);
    indices = round(linspace(1, n_points, min(5, n_points)));

    for idx = indices
        fprintf('\n--- 频率点: %.3f GHz ---\n', results.frequency(idx)/1e9);
        fprintf('S11 = %.2f dB ∠%.1f°\n', ...
                20*log10(abs(results.S11(idx))), angle(results.S11(idx))*180/pi);
        fprintf('S21 = %.2f dB ∠%.1f°\n', ...
                20*log10(abs(results.S21(idx))), angle(results.S21(idx))*180/pi);
        fprintf('\n反射系数 Γ = %.4f ∠%.1f°\n', ...
                abs(results.Gamma(idx)), angle(results.Gamma(idx))*180/pi);
        fprintf('传输系数 T = %.4f ∠%.1f°\n', ...
                abs(results.T(idx)), angle(results.T(idx))*180/pi);
        fprintf('\n相对磁导率:\n');
        fprintf('  μr = %.4f + j%.4f\n', ...
                real(results.mu_r(idx)), imag(results.mu_r(idx)));
        fprintf('  |μr| = %.4f\n', abs(results.mu_r(idx)));
        fprintf('\n相对介电常数:\n');
        fprintf('  εr = %.4f + j%.4f\n', ...
                real(results.epsilon_r(idx)), imag(results.epsilon_r(idx)));
        fprintf('  |εr| = %.4f\n', abs(results.epsilon_r(idx)));
        fprintf('\n损耗角正切: tanδ = %.6f\n', results.tan_delta(idx));
    end

    %% ========== 输出文件说明 ==========
    fprintf('\n\n========== 输出文件 ==========\n');
    fprintf('✓ Excel结果表格: results_from_s2p.xlsx\n');
    fprintf('✓ 结果图像: s2p_results_plot.png\n');
    fprintf('\nExcel表格包含:\n');
    fprintf('  - 所有频率点的完整数据\n');
    fprintf('  - S参数（dB格式）\n');
    fprintf('  - 反射/传输系数（幅度和相位）\n');
    fprintf('  - μr（实部、虚部、幅度）\n');
    fprintf('  - εr（实部、虚部、幅度）\n');
    fprintf('  - 损耗角正切\n');

    fprintf('\n图像包含6个子图:\n');
    fprintf('  1. S参数 vs 频率\n');
    fprintf('  2. 反射系数和传输系数 vs 频率\n');
    fprintf('  3. μr实部 vs 频率\n');
    fprintf('  4. μr虚部 vs 频率\n');
    fprintf('  5. εr实部 vs 频率\n');
    fprintf('  6. tanδ vs 频率\n');

    %% ========== 简单分析 ==========
    fprintf('\n\n========== 材料特性简析 ==========\n');

    % 判断磁性
    avg_mu_r = mean(real(results.mu_r));
    if abs(avg_mu_r - 1) < 0.1
        fprintf('✓ 材料类型: 非磁性材料 (μr ≈ 1)\n');
    else
        fprintf('✓ 材料类型: 磁性材料 (μr ≠ 1)\n');
    end

    % 介电常数
    avg_epsilon_r = mean(real(results.epsilon_r));
    fprintf('✓ 平均相对介电常数: εr'' ≈ %.2f\n', avg_epsilon_r);

    % 损耗评估
    avg_tan_delta = mean(results.tan_delta);
    fprintf('✓ 平均损耗角正切: tanδ ≈ %.4f\n', avg_tan_delta);

    if avg_tan_delta < 0.001
        fprintf('  → 低损耗材料\n');
    elseif avg_tan_delta < 0.01
        fprintf('  → 中等损耗材料\n');
    else
        fprintf('  → 高损耗材料\n');
    end

catch ME
    fprintf('\n❌ 错误: %s\n\n', ME.message);
    fprintf('请检查:\n');
    fprintf('  1. S2P文件路径是否正确\n');
    fprintf('  2. S2P文件格式是否标准\n');
    fprintf('  3. 样品厚度是否正确设置（单位: 米）\n');
    fprintf('  4. 频率范围是否在X波段内（8-12 GHz）\n\n');

    fprintf('常见问题:\n');
    fprintf('  • 如果提示"无法打开文件"，请检查文件路径\n');
    fprintf('  • 如果提示"未能读取到有效数据"，请检查S2P文件格式\n');
    fprintf('  • 如果计算结果异常，可能是样品厚度设置错误\n\n');

    rethrow(ME);
end

fprintf('\n\n========== 计算完成 ==========\n\n');

%% ========== 快速使用提示 ==========
fprintf('╔════════════════════════════════════════════════════╗\n');
fprintf('║                 快速使用方法                       ║\n');
fprintf('╠════════════════════════════════════════════════════╣\n');
fprintf('║  1. 修改 s2p_file 为你的S2P文件路径              ║\n');
fprintf('║  2. 修改 d 为你的样品厚度（单位: 米）            ║\n');
fprintf('║  3. 运行本脚本或直接调用:                         ║\n');
fprintf('║     results = calculate_from_s2p(file, d)         ║\n');
fprintf('╚════════════════════════════════════════════════════╝\n');

%% ========== 厚度参考 ==========
fprintf('\n常用厚度参考:\n');
fprintf('  0.5 mm  →  d = 0.5e-3\n');
fprintf('  1.0 mm  →  d = 1.0e-3\n');
fprintf('  2.0 mm  →  d = 2.0e-3\n');
fprintf('  3.0 mm  →  d = 3.0e-3\n');
fprintf('  5.0 mm  →  d = 5.0e-3\n');
fprintf('  10 mm   →  d = 10e-3\n');
