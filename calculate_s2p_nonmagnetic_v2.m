function [results] = calculate_s2p_nonmagnetic_v2(s2p_file, d, options)
% 非磁性材料专用计算V2（更智能的约束处理）
%
% 改进：
%   - 两阶段评分（先找可行解，再优化）
%   - 更好的回退策略
%   - 对问题频率点的特殊处理

    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;
    end

    if ~isfield(options, 'k_range')
        options.k_range = 20;  % 默认扩大到±20
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    if options.verbose
        fprintf('========== 非磁性材料参数计算 V2 ==========\n\n');
        fprintf('改进: 智能约束处理，更好的回退策略\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  样品厚度 d = %.2f mm\n', d * 1000);
        fprintf('  相位搜索范围: ±%d\n\n', options.k_range);
    end

    %% 读取S2P
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);
    lambda_0 = c ./ freq;

    %% 初始化
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    mu_r = ones(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);
    warnings_list = cell(n_points, 1);
    selected_k = zeros(n_points, 1);
    quality_score = zeros(n_points, 1);

    fprintf('开始计算...\n');

    for i = 1:n_points
        try
            %% 步骤1: 计算Gamma候选
            X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
            Gamma_p = X + sqrt(X^2 - 1);
            Gamma_m = X - sqrt(X^2 - 1);

            gamma_candidates = [Gamma_p; Gamma_m];

            %% 步骤2: 对每个Gamma和相位分支组合评分
            all_results = [];  % 存储所有候选结果

            for g_idx = 1:length(gamma_candidates)
                Gamma_test = gamma_candidates(g_idx);

                % 计算T
                T_test = (S11(i) + S21(i) - Gamma_test) / ...
                         (1 - (S11(i) + S21(i)) * Gamma_test);

                %% 对每个相位分支
                for k = -options.k_range:options.k_range
                    lgT = log(1/T_test) + 2*pi*1j*k;
                    inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
                    epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

                    %% 两阶段评分
                    % 阶段1: 硬约束（必须满足）
                    hard_score = 0;
                    hard_violations = 0;

                    % H1: εr'必须>0（放宽到0）
                    if real(epsilon_r_test) <= 0
                        hard_violations = hard_violations + 100;
                    else
                        hard_score = hard_score + 50;
                    end

                    % H2: |Γ|尽量≤1.1（放宽）
                    if abs(Gamma_test) > 1.2
                        hard_violations = hard_violations + 50;
                    elseif abs(Gamma_test) <= 1.0
                        hard_score = hard_score + 30;
                    end

                    % H3: |T|尽量≤1.1（放宽）
                    if abs(T_test) > 1.3
                        hard_violations = hard_violations + 100;
                    elseif abs(T_test) <= 1.0
                        hard_score = hard_score + 40;
                    else
                        hard_score = hard_score + 20 - (abs(T_test) - 1) * 50;
                    end

                    % 阶段2: 软约束（优化目标）
                    soft_score = 0;

                    % S1: εr'合理范围（1-100）
                    eps_real = real(epsilon_r_test);
                    if eps_real >= 1 && eps_real <= 50
                        soft_score = soft_score + 30;
                    elseif eps_real > 50 && eps_real <= 100
                        soft_score = soft_score + 10;
                    else
                        soft_score = soft_score - abs(eps_real - 25);
                    end

                    % S2: εr''>=0（损耗为正）
                    eps_imag = imag(epsilon_r_test);
                    if eps_imag >= 0
                        soft_score = soft_score + 25;
                    else
                        soft_score = soft_score - abs(eps_imag) * 30;
                    end

                    % S3: 连续性
                    if i > 1 && ~isnan(epsilon_r(i-1))
                        diff = abs(eps_real - real(epsilon_r(i-1)));
                        if diff < 2
                            soft_score = soft_score + 40;
                        elseif diff < 5
                            soft_score = soft_score + 20;
                        elseif diff < 10
                            soft_score = soft_score + 5;
                        else
                            soft_score = soft_score - diff * 2;
                        end
                    end

                    % S4: tanδ合理范围
                    if eps_real > 0.1
                        tan_d = abs(eps_imag / eps_real);
                        if tan_d >= 0.001 && tan_d <= 0.5
                            soft_score = soft_score + 15;
                        elseif tan_d > 0.5 && tan_d < 1.0
                            soft_score = soft_score + 5;
                        end
                    end

                    % 总分 = 硬约束分 - 违规惩罚 + 软约束分
                    total_score = hard_score - hard_violations + soft_score;

                    % 保存结果
                    candidate = struct();
                    candidate.Gamma = Gamma_test;
                    candidate.T = T_test;
                    candidate.epsilon_r = epsilon_r_test;
                    candidate.k = k;
                    candidate.g_idx = g_idx;
                    candidate.score = total_score;
                    candidate.hard_violations = hard_violations;

                    all_results = [all_results; candidate];
                end
            end

            %% 步骤3: 选择最佳结果（智能策略）
            if ~isempty(all_results)
                % 优先选择无硬约束违规的
                no_violation_idx = find([all_results.hard_violations] == 0);

                if ~isempty(no_violation_idx)
                    % 有无违规的解，从中选最高分
                    [~, best_idx] = max([all_results(no_violation_idx).score]);
                    best_result = all_results(no_violation_idx(best_idx));
                else
                    % 没有完美解，选择违规最少且分数最高的
                    [~, best_idx] = max([all_results.score]);
                    best_result = all_results(best_idx);
                    warnings_list{i} = sprintf('次优解(违规=%d)', best_result.hard_violations);
                end

                % 应用最佳结果
                Gamma(i) = best_result.Gamma;
                T(i) = best_result.T;
                epsilon_r(i) = best_result.epsilon_r;
                selected_k(i) = best_result.k;
                quality_score(i) = best_result.score;

                % 检查结果质量
                if abs(T(i)) > 1.05
                    if isempty(warnings_list{i})
                        warnings_list{i} = sprintf('|T|=%.3f', abs(T(i)));
                    else
                        warnings_list{i} = sprintf('%s; |T|=%.3f', warnings_list{i}, abs(T(i)));
                    end
                end

                if abs(Gamma(i)) > 1.05
                    if isempty(warnings_list{i})
                        warnings_list{i} = sprintf('|Γ|=%.3f', abs(Gamma(i)));
                    else
                        warnings_list{i} = sprintf('%s; |Γ|=%.3f', warnings_list{i}, abs(Gamma(i)));
                    end
                end

                if real(epsilon_r(i)) < 0.5
                    if isempty(warnings_list{i})
                        warnings_list{i} = sprintf('εr''=%.2f<1', real(epsilon_r(i)));
                    else
                        warnings_list{i} = sprintf('%s; εr''=%.2f', warnings_list{i}, real(epsilon_r(i)));
                    end
                end
            else
                % 极端情况：完全没有候选
                warnings_list{i} = '无任何候选解';
                epsilon_r(i) = NaN;
            end

            % 计算tanδ
            if abs(real(epsilon_r(i))) > 1e-6
                tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
            else
                tan_delta(i) = NaN;
            end

        catch ME
            warnings_list{i} = sprintf('错误: %s', ME.message);
            epsilon_r(i) = NaN;
        end

        % 进度
        if options.verbose && mod(i, round(n_points/10)) == 0
            fprintf('  进度: %d/%d (%.1f%%)\n', i, n_points, i/n_points*100);
        end
    end

    fprintf('计算完成！\n\n');

    %% 诊断报告
    if options.verbose
        fprintf('========== 诊断报告 ==========\n\n');

        n_warnings = sum(~cellfun(@isempty, warnings_list));
        fprintf('有警告的点: %d / %d (%.1f%%)\n\n', n_warnings, n_points, n_warnings/n_points*100);

        % 物理约束统计
        valid_Gamma = sum(abs(Gamma) <= 1.05);
        valid_T = sum(abs(T) <= 1.05);
        valid_eps = sum(real(epsilon_r) >= 0.5 & real(epsilon_r) < 100);
        valid_eps_imag = sum(imag(epsilon_r) >= -0.1);  % 允许小负值

        fprintf('物理约束符合率（放宽）:\n');
        fprintf('  |Γ| ≤ 1.05: %d/%d (%.1f%%)\n', valid_Gamma, n_points, valid_Gamma/n_points*100);
        fprintf('  |T| ≤ 1.05: %d/%d (%.1f%%)\n', valid_T, n_points, valid_T/n_points*100);
        fprintf('  εr'' ≥ 0.5: %d/%d (%.1f%%)\n', valid_eps, n_points, valid_eps/n_points*100);
        fprintf('  εr'''' ≥ -0.1: %d/%d (%.1f%%)\n', valid_eps_imag, n_points, valid_eps_imag/n_points*100);

        % 质量评分分布
        fprintf('\n质量评分分布:\n');
        high_quality = sum(quality_score > 100);
        medium_quality = sum(quality_score > 0 & quality_score <= 100);
        low_quality = sum(quality_score <= 0);

        fprintf('  高质量 (>100分): %d (%.1f%%)\n', high_quality, high_quality/n_points*100);
        fprintf('  中等质量 (0-100分): %d (%.1f%%)\n', medium_quality, medium_quality/n_points*100);
        fprintf('  低质量 (<0分): %d (%.1f%%)\n', low_quality, low_quality/n_points*100);

        % 推荐频率范围
        good_idx = (abs(T) <= 1.05) & (real(epsilon_r) >= 1) & (imag(epsilon_r) >= 0);
        if sum(good_idx) > 0
            fprintf('\n严格高质量点: %d/%d (%.1f%%)\n', sum(good_idx), n_points, sum(good_idx)/n_points*100);
            fprintf('推荐频率范围: %.3f - %.3f GHz\n', ...
                    min(freq(good_idx))/1e9, max(freq(good_idx))/1e9);

            fprintf('\n材料参数（高质量点）:\n');
            fprintf('  εr'' = %.3f ± %.3f\n', mean(real(epsilon_r(good_idx))), std(real(epsilon_r(good_idx))));
            fprintf('  εr'''' = %.3f ± %.3f\n', mean(imag(epsilon_r(good_idx))), std(imag(epsilon_r(good_idx))));
            fprintf('  tanδ = %.4f ± %.4f\n', mean(tan_delta(good_idx)), std(tan_delta(good_idx)));
        end
    end

    %% 保存结果
    S11_dB = 20*log10(abs(S11));
    S21_dB = 20*log10(abs(S21));

    output_file = 'results_nonmagnetic_v2.xlsx';
    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, selected_k, quality_score, warnings_list, ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_mag', ...
                                          'tan_delta', 'Selected_k', 'Quality_Score', 'Warnings'});

    writetable(output_table, output_file);
    fprintf('\n结果已保存到: %s\n', output_file);

    %% 绘图
    plot_results_v2(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, quality_score, good_idx);

    %% 返回
    results.frequency = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.mu_r = mu_r;
    results.epsilon_r = epsilon_r;
    results.tan_delta = tan_delta;
    results.selected_k = selected_k;
    results.quality_score = quality_score;
    results.warnings = warnings_list;
    results.good_idx = good_idx;

    fprintf('\n========== 完成 ==========\n');
end

function plot_results_v2(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, quality_score, good_idx)
    figure('Position', [100, 100, 1600, 1000], 'Name', '非磁性材料分析 V2');

    % 子图1: S参数
    subplot(3, 3, 1);
    plot(freq/1e9, S11_dB, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, S21_dB, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    legend('S_{11}', 'S_{21}');
    title('S参数');

    % 子图2: |Γ|和|T|
    subplot(3, 3, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1.0, 'k--', '理想上限', 'LineWidth', 1);
    yline(1.05, 'k:', '放宽上限', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|', 'Location', 'best');
    title('反射/传输系数');

    % 子图3: εr实部
    subplot(3, 3, 3);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 1.5);
    hold on;
    yline(1.0, 'r--', '真空下限', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, real(epsilon_r(good_idx)), 'g.', 'MarkerSize', 8);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部）');

    % 子图4: εr虚部
    subplot(3, 3, 4);
    plot(freq/1e9, imag(epsilon_r), 'r-', 'LineWidth', 1.5);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, imag(epsilon_r(good_idx)), 'g.', 'MarkerSize', 8);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''''');
    title('相对介电常数（虚部）');

    % 子图5: tanδ
    subplot(3, 3, 5);
    valid_tan = abs(tan_delta) < 10;
    plot(freq(valid_tan)/1e9, tan_delta(valid_tan), 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    % 子图6: 质量评分
    subplot(3, 3, 6);
    plot(freq/1e9, quality_score, 'k-', 'LineWidth', 2);
    hold on;
    yline(100, 'g--', '高质量阈值', 'LineWidth', 1);
    yline(0, 'r--', '及格线', 'LineWidth', 1);
    fill([freq'/1e9; flipud(freq'/1e9)], ...
         [max(quality_score, 0); zeros(size(quality_score))], ...
         'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('质量评分');
    title('数据质量评分');

    % 子图7: |Γ|详细
    subplot(3, 3, 7);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    violation_idx = abs(Gamma) > 1.05;
    if sum(violation_idx) > 0
        plot(freq(violation_idx)/1e9, abs(Gamma(violation_idx)), 'ro', 'MarkerSize', 6);
    end
    yline(1.0, 'k--', '理想', 'LineWidth', 1);
    yline(1.05, 'r--', '放宽', 'LineWidth', 1);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|\Gamma|');
    title('反射系数详细（红圈=超限）');
    legend('|\Gamma|', '超限点', 'Location', 'best');

    % 子图8: |T|详细
    subplot(3, 3, 8);
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    hold on;
    violation_idx = abs(T) > 1.05;
    if sum(violation_idx) > 0
        plot(freq(violation_idx)/1e9, abs(T(violation_idx)), 'ro', 'MarkerSize', 6);
    end
    yline(1.0, 'k--', '理想', 'LineWidth', 1);
    yline(1.05, 'r--', '放宽', 'LineWidth', 1);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|T|');
    title('传输系数详细（红圈=超限）');
    legend('|T|', '超限点', 'Location', 'best');

    % 子图9: 问题点分布
    subplot(3, 3, 9);
    problem_flags = zeros(size(freq));
    problem_flags(abs(Gamma) > 1.05) = problem_flags(abs(Gamma) > 1.05) + 1;
    problem_flags(abs(T) > 1.05) = problem_flags(abs(T) > 1.05) + 1;
    problem_flags(real(epsilon_r) < 1) = problem_flags(real(epsilon_r) < 1) + 1;
    problem_flags(imag(epsilon_r) < -0.1) = problem_flags(imag(epsilon_r) < -0.1) + 1;

    bar(freq/1e9, problem_flags, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'none');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('违规约束数量');
    title('问题点分布 (0=完美)');
    ylim([0, 4.5]);

    saveas(gcf, 'results_nonmagnetic_v2_plot.png');
    fprintf('图像已保存到: results_nonmagnetic_v2_plot.png\n');
end
