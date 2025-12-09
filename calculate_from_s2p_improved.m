function [results] = calculate_from_s2p_improved(s2p_file, d, options)
% 改进版S2P材料参数计算（修复相位解缠和数值稳定性问题）
%
% 输入参数:
%   s2p_file: S2P文件路径
%   d: 样品厚度 (单位: 米)
%   options: 可选参数结构体
%       - lambda_c: 截止波长 (m)，默认 45.7e-3 (X波段)
%       - max_iterations: 相位解缠最大迭代次数，默认 10
%       - verbose: 是否显示详细信息，默认 true
%
% 改进内容:
%   1. 改进的相位解缠算法
%   2. 多次迭代选择最优Γ符号
%   3. 物理约束检查和修正
%   4. 异常点诊断

    %% ========== 参数设置 ==========
    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;  % X波段默认
    end

    if ~isfield(options, 'max_iterations')
        options.max_iterations = 10;
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    % 物理常数
    c = 299792458;
    mu_0 = 4*pi*1e-7;
    epsilon_0 = 8.854187817e-12;

    lambda_c = options.lambda_c;
    waveguide_a = lambda_c / 2;

    epsilon_i = 1;
    mu_i = 1;

    if options.verbose
        fprintf('========== 改进版材料参数计算 ==========\n\n');
        fprintf('固定参数:\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  波导宽边 a = %.2f mm\n', waveguide_a * 1000);
        fprintf('  样品厚度 d = %.2f mm\n', d * 1000);
        fprintf('\n');
    end

    %% ========== 读取S2P文件 ==========
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);

    % 检查频率范围
    freq_ghz_min = min(freq) / 1e9;
    freq_ghz_max = max(freq) / 1e9;

    if options.verbose
        fprintf('频率范围: %.3f - %.3f GHz (%d 个点)\n\n', ...
                freq_ghz_min, freq_ghz_max, n_points);
    end

    %% ========== 改进的计算算法 ==========
    lambda_0 = c ./ freq;

    % 初始化结果数组
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    mu_r = zeros(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);
    lambda_g = zeros(n_points, 1);

    % 记录警告信息
    warnings = cell(n_points, 1);

    fprintf('开始计算（改进算法）...\n');

    for i = 1:n_points
        try
            % 方程 (2-40): 计算 X
            X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));

            % 方程 (2-41): 计算两个可能的 Gamma
            sqrt_term = sqrt(X^2 - 1);
            Gamma_p = X + sqrt_term;
            Gamma_m = X - sqrt_term;

            % 改进的符号选择策略
            % 1. 首先检查 |Gamma| <= 1 约束
            mag_p = abs(Gamma_p);
            mag_m = abs(Gamma_m);

            valid_p = (mag_p <= 1.01);  % 允许小误差
            valid_m = (mag_m <= 1.01);

            % 2. 如果两个都有效，选择幅度较小的
            if valid_p && valid_m
                if mag_m < mag_p
                    Gamma(i) = Gamma_m;
                else
                    Gamma(i) = Gamma_p;
                end
            elseif valid_p
                Gamma(i) = Gamma_p;
            elseif valid_m
                Gamma(i) = Gamma_m;
            else
                % 两个都不满足，选择更接近1的
                if mag_m < mag_p
                    Gamma(i) = Gamma_m / mag_m * 0.99;  % 归一化到 <1
                else
                    Gamma(i) = Gamma_p / mag_p * 0.99;
                end
                warnings{i} = 'Gamma超出约束，已修正';
            end

            % 计算传输系数 T
            T(i) = (S11(i) + S21(i) - Gamma(i)) / (1 - (S11(i) + S21(i)) * Gamma(i));

            % 检查 |T| 约束
            if abs(T(i)) > 1.01
                if isempty(warnings{i})
                    warnings{i} = sprintf('|T|=%.3f>1', abs(T(i)));
                else
                    warnings{i} = sprintf('%s; |T|=%.3f>1', warnings{i}, abs(T(i)));
                end

                % 尝试另一个Gamma
                if valid_p && ~valid_m
                    Gamma_alt = Gamma_m;
                else
                    Gamma_alt = Gamma_p;
                end
                T_alt = (S11(i) + S21(i) - Gamma_alt) / (1 - (S11(i) + S21(i)) * Gamma_alt);

                if abs(T_alt) < abs(T(i))
                    Gamma(i) = Gamma_alt;
                    T(i) = T_alt;
                    warnings{i} = sprintf('%s (已切换Gamma)', warnings{i});
                end
            end

            % 改进的相位解缠
            % 计算 ln(1/T)，需要处理多值性
            log_inv_T = -log(T(i));  % ln(1/T) = -ln(T)

            % 方法：尝试不同的相位分支，选择物理合理的结果
            best_branch = 0;
            best_score = -inf;

            for branch = -options.max_iterations:options.max_iterations
                % 调整相位
                log_test = log_inv_T + 2*pi*1j*branch;

                % 计算对应的 1/lambda^2
                inv_lambda_sq_test = -(1/(2*pi*d) * log_test)^2;

                % 计算 mu_r 和 epsilon_r
                term1 = sqrt(epsilon_i * mu_i / lambda_0(i)^2 - 1/lambda_c^2);
                Lambda = d;

                mu_r_test = (1 + Gamma(i)) * mu_i / (Lambda * (1 - Gamma(i)) * term1);

                term2 = (inv_lambda_sq_test + 1/lambda_c^2) * lambda_0(i)^2;
                epsilon_r_test = term2 / mu_r_test;

                % 评分标准（物理合理性）
                score = 0;

                % 1. epsilon_r 实部应为正
                if real(epsilon_r_test) > 0.5
                    score = score + 10;
                end

                % 2. mu_r 实部应接近1（对非磁性材料）
                if abs(real(mu_r_test) - 1) < 1
                    score = score + 5;
                end

                % 3. epsilon_r 虚部应为正（损耗）
                if imag(epsilon_r_test) > 0
                    score = score + 5;
                end

                % 4. 值应该合理
                if real(epsilon_r_test) < 100 && abs(real(mu_r_test)) < 10
                    score = score + 5;
                end

                % 5. 连续性（如果不是第一个点）
                if i > 1
                    continuity_score = -abs(real(epsilon_r_test) - real(epsilon_r(i-1)));
                    score = score + continuity_score * 0.1;
                end

                if score > best_score
                    best_score = score;
                    best_branch = branch;
                    mu_r(i) = mu_r_test;
                    epsilon_r(i) = epsilon_r_test;
                    lambda_g(i) = 1 / sqrt(inv_lambda_sq_test);
                end
            end

            % 计算损耗角正切
            if abs(real(epsilon_r(i))) > 1e-6
                tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
            else
                tan_delta(i) = NaN;
                if isempty(warnings{i})
                    warnings{i} = 'εr实部接近0';
                else
                    warnings{i} = sprintf('%s; εr实部接近0', warnings{i});
                end
            end

            % 最终物理检查
            if real(epsilon_r(i)) < 0
                if isempty(warnings{i})
                    warnings{i} = sprintf('εr实部<0 (分支%d)', best_branch);
                else
                    warnings{i} = sprintf('%s; εr实部<0 (分支%d)', warnings{i}, best_branch);
                end
            end

            if abs(real(mu_r(i))) > 20
                if isempty(warnings{i})
                    warnings{i} = 'μr异常大';
                else
                    warnings{i} = sprintf('%s; μr异常大', warnings{i});
                end
            end

        catch ME
            warnings{i} = sprintf('计算错误: %s', ME.message);
            mu_r(i) = NaN;
            epsilon_r(i) = NaN;
            tan_delta(i) = NaN;
        end

        % 进度显示
        if options.verbose && mod(i, round(n_points/10)) == 0
            fprintf('  进度: %d/%d (%.1f%%)\n', i, n_points, i/n_points*100);
        end
    end

    fprintf('计算完成！\n\n');

    %% ========== 诊断和报告 ==========
    if options.verbose
        fprintf('========== 诊断报告 ==========\n\n');

        % 统计警告
        n_warnings = sum(~cellfun(@isempty, warnings));
        if n_warnings > 0
            fprintf('⚠️  发现 %d 个问题点:\n', n_warnings);
            for i = 1:n_points
                if ~isempty(warnings{i})
                    fprintf('  频率 %.3f GHz: %s\n', freq(i)/1e9, warnings{i});
                end
            end
            fprintf('\n');
        end

        % 物理约束检查
        fprintf('物理约束检查:\n');

        violation_gamma = sum(abs(Gamma) > 1.01);
        violation_T = sum(abs(T) > 1.01);
        violation_epsilon = sum(real(epsilon_r) < 0);
        violation_mu = sum(abs(real(mu_r)) > 10 | real(mu_r) < 0);

        fprintf('  |Γ| > 1: %d 个点\n', violation_gamma);
        fprintf('  |T| > 1: %d 个点\n', violation_T);
        fprintf('  εr'' < 0: %d 个点\n', violation_epsilon);
        fprintf('  μr'' 异常: %d 个点\n\n', violation_mu);

        % 统计信息
        fprintf('========== 计算结果统计 ==========\n\n');

        % 排除NaN和异常值
        valid_idx = ~isnan(real(mu_r)) & abs(real(mu_r)) < 10;

        if sum(valid_idx) > 0
            fprintf('相对磁导率 μr (有效点 %d/%d):\n', sum(valid_idx), n_points);
            fprintf('  实部范围: %.4f ~ %.4f\n', min(real(mu_r(valid_idx))), max(real(mu_r(valid_idx))));
            fprintf('  平均值: %.4f + j%.4f\n', mean(real(mu_r(valid_idx))), mean(imag(mu_r(valid_idx))));
        end

        valid_idx = ~isnan(real(epsilon_r)) & real(epsilon_r) > 0 & real(epsilon_r) < 100;

        if sum(valid_idx) > 0
            fprintf('\n相对介电常数 εr (有效点 %d/%d):\n', sum(valid_idx), n_points);
            fprintf('  实部范围: %.4f ~ %.4f\n', min(real(epsilon_r(valid_idx))), max(real(epsilon_r(valid_idx))));
            fprintf('  平均值: %.4f + j%.4f\n', mean(real(epsilon_r(valid_idx))), mean(imag(epsilon_r(valid_idx))));
        end

        valid_idx = ~isnan(tan_delta) & abs(tan_delta) < 10;

        if sum(valid_idx) > 0
            fprintf('\n损耗角正切 tanδ (有效点 %d/%d):\n', sum(valid_idx), n_points);
            fprintf('  范围: %.6f ~ %.6f\n', min(tan_delta(valid_idx)), max(tan_delta(valid_idx)));
            fprintf('  平均值: %.6f\n', mean(tan_delta(valid_idx)));
        end
    end

    %% ========== 保存结果 ==========
    S11_dB = 20 * log10(abs(S11));
    S21_dB = 20 * log10(abs(S21));

    output_file = 'results_improved.xlsx';

    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(mu_r), imag(mu_r), abs(mu_r), ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, warnings, ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'mu_r_real', 'mu_r_imag', 'mu_r_mag', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_mag', ...
                                          'tan_delta', 'Warnings'});

    writetable(output_table, output_file);

    if options.verbose
        fprintf('\n结果已保存到: %s\n', output_file);
    end

    %% ========== 绘图 ==========
    plot_improved_results(freq, S11_dB, S21_dB, Gamma, T, mu_r, epsilon_r, tan_delta);

    %% ========== 返回结果 ==========
    results.frequency = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.mu_r = mu_r;
    results.epsilon_r = epsilon_r;
    results.tan_delta = tan_delta;
    results.lambda_0 = lambda_0;
    results.lambda_c = lambda_c;
    results.lambda_g = lambda_g;
    results.warnings = warnings;

    fprintf('\n========== 完成 ==========\n');
end

function plot_improved_results(freq, S11_dB, S21_dB, Gamma, T, mu_r, epsilon_r, tan_delta)
    % 绘制改进后的结果

    figure('Position', [100, 100, 1400, 900], 'Name', '改进版材料参数分析结果');

    % 子图1: S参数
    subplot(3, 2, 1);
    plot(freq/1e9, S11_dB, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, S21_dB, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    legend('S_{11}', 'S_{21}');
    title('S参数');

    % 子图2: |Γ| 和 |T| (带约束线)
    subplot(3, 2, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1.0, 'k--', '物理上限', 'LineWidth', 1);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|', 'Location', 'best');
    title('反射/传输系数（应 ≤ 1）');
    ylim([0, max(1.5, max([abs(Gamma); abs(T)]))]);

    % 子图3: μr 实部
    subplot(3, 2, 3);
    plot(freq/1e9, real(mu_r), 'b-', 'LineWidth', 2);
    hold on;
    yline(1.0, 'k--', '非磁性材料', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\mu_r''');
    title('相对磁导率（实部）');

    % 子图4: εr 实部
    subplot(3, 2, 4);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 2);
    hold on;
    yline(0, 'r--', '物理下限', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部，应 > 0）');

    % 子图5: εr 虚部
    subplot(3, 2, 5);
    plot(freq/1e9, imag(epsilon_r), 'r-', 'LineWidth', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''''');
    title('相对介电常数（虚部，损耗）');

    % 子图6: tanδ
    subplot(3, 2, 6);
    % 只绘制合理范围内的值
    valid_idx = abs(tan_delta) < 10;
    plot(freq(valid_idx)/1e9, tan_delta(valid_idx), 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    saveas(gcf, 'results_improved_plot.png');
    fprintf('图像已保存到: results_improved_plot.png\n');
end
