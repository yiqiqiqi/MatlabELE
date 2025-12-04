%% 改进版S2P材料参数计算 - 使用示例
%
% 本脚本演示如何使用改进版的 calculate_from_s2p_improved 函数
% 解决了原版本的相位解缠和数值稳定性问题
%
% 改进内容:
%   1. 智能相位解缠（尝试多个分支，选择最合理的结果）
%   2. 改进的Γ符号选择（考虑|T|约束）
%   3. 物理约束检查和自动修正
%   4. 详细的诊断报告

clear all;
close all;
clc;

%% ========== 参数设置 ==========

% S2P文件路径
s2p_file = 'your_material_test.s2p';  % ← 修改为你的文件路径

% 样品厚度 (单位: 米)
d = 1.4e-3;  % 1.4mm

% 可选参数
options.lambda_c = 45.7e-3;  % 截止波长 (m)，X波段默认
options.max_iterations = 10;  % 相位解缠最大尝试次数
options.verbose = true;  % 显示详细信息

%% ========== 运行改进版计算 ==========

fprintf('╔════════════════════════════════════════════════════╗\n');
fprintf('║        改进版材料参数计算（修复版）              ║\n');
fprintf('╚════════════════════════════════════════════════════╝\n\n');

try
    % 调用改进版函数
    results = calculate_from_s2p_improved(s2p_file, d, options);

    %% ========== 结果分析 ==========
    fprintf('\n╔════════════════════════════════════════════════════╗\n');
    fprintf('║                结果质量评估                       ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');

    % 找出有效数据点
    valid_epsilon = real(results.epsilon_r) > 0 & real(results.epsilon_r) < 100;
    valid_mu = abs(real(results.mu_r)) < 10 & real(results.mu_r) > 0;
    valid_T = abs(results.T) <= 1.01;
    valid_Gamma = abs(results.Gamma) <= 1.01;

    fprintf('数据质量:\n');
    fprintf('  总数据点: %d\n', length(results.frequency));
    fprintf('  |Γ|≤1 的点: %d (%.1f%%)\n', sum(valid_Gamma), sum(valid_Gamma)/length(results.frequency)*100);
    fprintf('  |T|≤1 的点: %d (%.1f%%)\n', sum(valid_T), sum(valid_T)/length(results.frequency)*100);
    fprintf('  εr合理的点: %d (%.1f%%)\n', sum(valid_epsilon), sum(valid_epsilon)/length(results.frequency)*100);
    fprintf('  μr合理的点: %d (%.1f%%)\n', sum(valid_mu), sum(valid_mu)/length(results.frequency)*100);

    % 选择"好"的频率范围
    good_points = valid_epsilon & valid_mu & valid_T & valid_Gamma;

    if sum(good_points) > 0
        fprintf('\n✓ 找到 %d 个高质量数据点\n', sum(good_points));
        fprintf('\n推荐使用的频率范围:\n');
        good_freq = results.frequency(good_points);
        fprintf('  %.3f - %.3f GHz\n', min(good_freq)/1e9, max(good_freq)/1e9);

        fprintf('\n在此频率范围内的材料参数:\n');
        fprintf('  相对磁导率 μr: %.3f + j%.3f\n', ...
                mean(real(results.mu_r(good_points))), ...
                mean(imag(results.mu_r(good_points))));
        fprintf('  相对介电常数 εr: %.3f + j%.3f\n', ...
                mean(real(results.epsilon_r(good_points))), ...
                mean(imag(results.epsilon_r(good_points))));
        fprintf('  损耗角正切 tanδ: %.4f\n', ...
                mean(results.tan_delta(good_points)));

        % 材料类型判断
        fprintf('\n材料特性分析:\n');
        avg_mu = mean(real(results.mu_r(good_points)));
        if abs(avg_mu - 1) < 0.2
            fprintf('  ✓ 非磁性材料 (μr ≈ 1)\n');
        else
            fprintf('  ✓ 磁性材料 (μr ≠ 1)\n');
        end

        avg_epsilon = mean(real(results.epsilon_r(good_points)));
        fprintf('  ✓ 相对介电常数约为 %.2f\n', avg_epsilon);

        avg_tan_delta = mean(results.tan_delta(good_points));
        if avg_tan_delta < 0.01
            fprintf('  ✓ 低损耗材料 (tanδ < 0.01)\n');
        elseif avg_tan_delta < 0.1
            fprintf('  ✓ 中等损耗材料\n');
        else
            fprintf('  ✓ 高损耗材料\n');
        end
    else
        fprintf('\n⚠️  警告：未找到高质量数据点\n');
        fprintf('建议检查:\n');
        fprintf('  1. S2P文件是否正确\n');
        fprintf('  2. 样品厚度是否准确\n');
        fprintf('  3. 波导参数是否匹配\n');
        fprintf('  4. 测量校准是否良好\n');
    end

    %% ========== 输出文件说明 ==========
    fprintf('\n\n╔════════════════════════════════════════════════════╗\n');
    fprintf('║                输出文件                           ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');
    fprintf('✓ Excel结果表格: results_improved.xlsx\n');
    fprintf('  - 包含所有计算结果\n');
    fprintf('  - 包含每个点的警告信息\n');
    fprintf('  - 可用于进一步分析\n\n');
    fprintf('✓ 结果图像: results_improved_plot.png\n');
    fprintf('  - 6个子图展示所有关键参数\n');
    fprintf('  - 标注了物理约束线\n');
    fprintf('  - 方便快速评估结果\n');

    %% ========== 与原版对比 ==========
    fprintf('\n\n╔════════════════════════════════════════════════════╗\n');
    fprintf('║            改进版 vs 原版对比                     ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');
    fprintf('改进版特点:\n');
    fprintf('  ✓ 自动尝试多个相位分支，选择最合理的结果\n');
    fprintf('  ✓ 智能Γ符号选择，考虑|T|≤1约束\n');
    fprintf('  ✓ 物理约束检查和自动修正\n');
    fprintf('  ✓ 详细的诊断报告\n');
    fprintf('  ✓ 标记问题点，方便排查\n\n');

    fprintf('建议:\n');
    fprintf('  • 优先使用"高质量数据点"范围内的结果\n');
    fprintf('  • 如果问题点较多，检查测量设置\n');
    fprintf('  • 可以尝试调整 options.max_iterations 增加搜索范围\n');

catch ME
    fprintf('\n❌ 错误: %s\n\n', ME.message);
    fprintf('可能的原因:\n');
    fprintf('  1. S2P文件路径不正确\n');
    fprintf('  2. S2P文件格式不标准\n');
    fprintf('  3. 样品厚度设置错误\n');
    fprintf('  4. 内存不足\n\n');
    rethrow(ME);
end

fprintf('\n\n========== 计算完成 ==========\n\n');

%% ========== 对比原版结果（如果存在）==========
if exist('results_from_s2p.xlsx', 'file')
    fprintf('\n提示: 检测到原版计算结果文件\n');
    fprintf('可以对比查看:\n');
    fprintf('  原版: results_from_s2p.xlsx\n');
    fprintf('  改进版: results_improved.xlsx\n');
end
