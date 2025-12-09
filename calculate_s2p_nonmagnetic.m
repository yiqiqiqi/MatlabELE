function [results] = calculate_s2p_nonmagnetic(s2p_file, d, options)
% 非磁性材料专用计算（μr固定为1）
% 简化版，专门针对塑料、陶瓷等非磁性材料
%
% 输入参数:
%   s2p_file: S2P文件路径
%   d: 样品厚度 (m)
%   options: 可选参数
%       - lambda_c: 截止波长 (m)，默认 45.7e-3
%       - k_range: 相位分支搜索范围，默认 15
%       - verbose: 显示详细信息，默认 true

    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;
    end

    if ~isfield(options, 'k_range')
        options.k_range = 15;
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    if options.verbose
        fprintf('========== 非磁性材料参数计算 (μr=1) ==========\n\n');
        fprintf('固定参数:\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  样品厚度 d = %.2f mm\n', d * 1000);
        fprintf('  假设: μr = 1 (非磁性材料)\n\n');
    end

    %% 读取S2P文件
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);
    lambda_0 = c ./ freq;

    %% 初始化
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    mu_r = ones(n_points, 1);  % 固定为1
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);
    warnings_list = cell(n_points, 1);
    selected_k = zeros(n_points, 1);

    fprintf('开始计算...\n');

    for i = 1:n_points
        try
            %% 步骤1: 计算Gamma
            X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
            Gamma_p = X + sqrt(X^2 - 1);
            Gamma_m = X - sqrt(X^2 - 1);

            mag_p = abs(Gamma_p);
            mag_m = abs(Gamma_m);

            % 选择策略：优先选|Gamma|小的
            gamma_candidates = [];
            if mag_p <= 1.05
                gamma_candidates = [gamma_candidates; Gamma_p];
            end
            if mag_m <= 1.05
                gamma_candidates = [gamma_candidates; Gamma_m];
            end

            if isempty(gamma_candidates)
                % 都不满足，选小的并归一化
                if mag_p < mag_m
                    Gamma(i) = Gamma_p / mag_p * 0.98;
                else
                    Gamma(i) = Gamma_m / mag_m * 0.98;
                end
                warnings_list{i} = '|Γ|>1已归一化';
            else
                % 选幅度最小的
                [~, idx] = min(abs(gamma_candidates));
                Gamma(i) = gamma_candidates(idx);
            end

            %% 步骤2: 对每个Gamma候选，尝试所有相位分支
            best_score = -inf;
            best_result = struct();

            % 尝试两个Gamma（如果都valid）
            gamma_to_test = Gamma(i);
            if length(gamma_candidates) > 1
                gamma_to_test = gamma_candidates;
            end

            for g_idx = 1:length(gamma_to_test)
                Gamma_test = gamma_to_test(g_idx);

                % 计算T
                T_test = (S11(i) + S21(i) - Gamma_test) / ...
                         (1 - (S11(i) + S21(i)) * Gamma_test);

                % 如果|T|>1.2，直接跳过这个Gamma
                if abs(T_test) > 1.2
                    continue;
                end

                %% 对每个相位分支k进行测试
                for k = -options.k_range:options.k_range
                    % 计算log(1/T) + 2πik
                    lgT = log(1/T_test) + 2*pi*1j*k;

                    % 计算1/λ²
                    inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;

                    % 计算εr（假设μr=1）
                    epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

                    %% 物理合理性评分
                    score = 0;

                    % 1. εr'必须>1（真空是下限）
                    eps_real = real(epsilon_r_test);
                    if eps_real < 1.0
                        continue;  % 直接跳过
                    end
                    score = score + 50;  % 基础分

                    % 2. εr'合理范围（2-50）
                    if eps_real >= 2 && eps_real <= 50
                        score = score + 30;
                    elseif eps_real < 2
                        score = score - (2 - eps_real) * 20;
                    elseif eps_real > 50
                        score = score - (eps_real - 50) * 2;
                    end

                    % 3. εr''必须>=0（损耗不能为负）
                    eps_imag = imag(epsilon_r_test);
                    if eps_imag < 0
                        score = score - abs(eps_imag) * 50;  // 强烈惩罚
                    else
                        score = score + 20;
                    end

                    % 4. |T|越接近1以下越好
                    T_mag = abs(T_test);
                    if T_mag <= 1.0
                        score = score + 30;
                    else
                        score = score - (T_mag - 1) * 100;  // 强烈惩罚
                    end

                    % 5. 连续性（与前一点接近）
                    if i > 1 && ~isnan(epsilon_r(i-1))
                        diff = abs(eps_real - real(epsilon_r(i-1)));
                        if diff < 1
                            score = score + 25;
                        elseif diff < 3
                            score = score + 10;
                        elseif diff < 5
                            score = score + 5;
                        else
                            score = score - diff * 3;
                        end

                        % 虚部连续性
                        diff_imag = abs(eps_imag - imag(epsilon_r(i-1)));
                        if diff_imag < 1
                            score = score + 10;
                        end
                    end

                    % 6. tanδ合理范围（0.001-0.5）
                    if eps_real > 0.1
                        tan_d = eps_imag / eps_real;
                        if tan_d >= 0.001 && tan_d <= 0.5
                            score = score + 15;
                        end
                    end

                    %% 记录最佳结果
                    if score > best_score
                        best_score = score;
                        best_result.Gamma = Gamma_test;
                        best_result.T = T_test;
                        best_result.epsilon_r = epsilon_r_test;
                        best_result.k = k;
                        best_result.score = score;
                    end
                end
            end

            %% 使用最佳结果
            if ~isempty(fieldnames(best_result))
                Gamma(i) = best_result.Gamma;
                T(i) = best_result.T;
                epsilon_r(i) = best_result.epsilon_r;
                selected_k(i) = best_result.k;

                % 检查最终结果
                if abs(T(i)) > 1.01
                    warnings_list{i} = sprintf('|T|=%.3f>1', abs(T(i)));
                end
                if real(epsilon_r(i)) < 1
                    if isempty(warnings_list{i})
                        warnings_list{i} = 'εr''<1';
                    else
                        warnings_list{i} = sprintf('%s; εr''<1', warnings_list{i});
                    end
                end
            else
                % 没找到合理解，使用默认
                Gamma(i) = Gamma_p;
                T(i) = (S11(i) + S21(i) - Gamma(i)) / ...
                       (1 - (S11(i) + S21(i)) * Gamma(i));
                lgT = log(1/T(i));
                inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
                epsilon_r(i) = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;
                warnings_list{i} = '无合理解，使用默认';
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
        fprintf('问题点数量: %d / %d\n\n', n_warnings, n_points);

        % 物理约束统计
        valid_Gamma = sum(abs(Gamma) <= 1.01);
        valid_T = sum(abs(T) <= 1.01);
        valid_eps = sum(real(epsilon_r) >= 1 & real(epsilon_r) < 100);
        valid_eps_imag = sum(imag(epsilon_r) >= 0);

        fprintf('物理约束符合率:\n');
        fprintf('  |Γ| ≤ 1: %d/%d (%.1f%%)\n', valid_Gamma, n_points, valid_Gamma/n_points*100);
        fprintf('  |T| ≤ 1: %d/%d (%.1f%%)\n', valid_T, n_points, valid_T/n_points*100);
        fprintf('  εr'' ≥ 1: %d/%d (%.1f%%)\n', valid_eps, n_points, valid_eps/n_points*100);
        fprintf('  εr'''' ≥ 0: %d/%d (%.1f%%)\n', valid_eps_imag, n_points, valid_eps_imag/n_points*100);

        % 统计信息
        good_idx = valid_T & valid_eps & valid_eps_imag;
        if sum(good_idx) > 0
            fprintf('\n高质量数据点: %d/%d (%.1f%%)\n', sum(good_idx), n_points, sum(good_idx)/n_points*100);
            fprintf('\n推荐使用频率范围: %.3f - %.3f GHz\n', ...
                    min(freq(good_idx))/1e9, max(freq(good_idx))/1e9);

            fprintf('\n材料参数（高质量点平均）:\n');
            fprintf('  εr'' = %.3f ± %.3f\n', mean(real(epsilon_r(good_idx))), std(real(epsilon_r(good_idx))));
            fprintf('  εr'''' = %.3f ± %.3f\n', mean(imag(epsilon_r(good_idx))), std(imag(epsilon_r(good_idx))));
            fprintf('  tanδ = %.4f ± %.4f\n', mean(tan_delta(good_idx)), std(tan_delta(good_idx)));
        end
    end

    %% 保存结果
    S11_dB = 20*log10(abs(S11));
    S21_dB = 20*log10(abs(S21));

    output_file = 'results_nonmagnetic.xlsx';
    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, selected_k, warnings_list, ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_mag', ...
                                          'tan_delta', 'Selected_k_branch', 'Warnings'});

    writetable(output_table, output_file);
    fprintf('\n结果已保存到: %s\n', output_file);

    %% 绘图
    plot_nonmagnetic_results(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, good_idx);

    %% 返回结果
    results.frequency = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.mu_r = mu_r;
    results.epsilon_r = epsilon_r;
    results.tan_delta = tan_delta;
    results.selected_k = selected_k;
    results.warnings = warnings_list;
    results.good_idx = good_idx;

    fprintf('\n========== 完成 ==========\n');
end

function plot_nonmagnetic_results(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, good_idx)
    figure('Position', [100, 100, 1600, 900], 'Name', '非磁性材料分析结果');

    % 子图1: S参数
    subplot(2, 3, 1);
    plot(freq/1e9, S11_dB, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, S21_dB, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    legend('S_{11}', 'S_{21}');
    title('S参数');

    % 子图2: |Γ|和|T|（带约束线）
    subplot(2, 3, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1.0, 'k--', '物理上限', 'LineWidth', 1);
    % 标记好点
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, abs(Gamma(good_idx)), 'b.', 'MarkerSize', 3);
        plot(freq(good_idx)/1e9, abs(T(good_idx)), 'r.', 'MarkerSize', 3);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|', 'Location', 'best');
    title('反射/传输系数 (应≤1)');
    ylim([0, min(1.5, max([abs(Gamma); abs(T)]))]);

    % 子图3: εr实部
    subplot(2, 3, 3);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 1.5);
    hold on;
    yline(1.0, 'r--', '物理下限 (真空)', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, real(epsilon_r(good_idx)), 'g.', 'MarkerSize', 5);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部，应≥1）');

    % 子图4: εr虚部
    subplot(2, 3, 4);
    plot(freq/1e9, imag(epsilon_r), 'r-', 'LineWidth', 1.5);
    hold on;
    yline(0, 'k--', '物理下限', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, imag(epsilon_r(good_idx)), 'g.', 'MarkerSize', 5);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''''');
    title('相对介电常数（虚部，应≥0）');

    % 子图5: tanδ
    subplot(2, 3, 5);
    valid_tan = abs(tan_delta) < 5;
    plot(freq(valid_tan)/1e9, tan_delta(valid_tan), 'k-', 'LineWidth', 1.5);
    hold on;
    if exist('good_idx', 'var')
        good_tan = good_idx & valid_tan;
        plot(freq(good_tan)/1e9, tan_delta(good_tan), 'g.', 'MarkerSize', 5);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    % 子图6: 数据质量指示
    subplot(2, 3, 6);
    quality = zeros(size(freq));
    quality(abs(Gamma) <= 1.01) = quality(abs(Gamma) <= 1.01) + 1;
    quality(abs(T) <= 1.01) = quality(abs(T) <= 1.01) + 1;
    quality(real(epsilon_r) >= 1) = quality(real(epsilon_r) >= 1) + 1;
    quality(imag(epsilon_r) >= 0) = quality(imag(epsilon_r) >= 0) + 1;

    plot(freq/1e9, quality, 'k-', 'LineWidth', 2);
    hold on;
    fill([freq'/1e9; flipud(freq'/1e9)], [quality; zeros(size(quality))], ...
         'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('通过约束数量');
    title('数据质量评分 (0-4分)');
    ylim([0, 4.5]);
    yticks(0:4);

    saveas(gcf, 'results_nonmagnetic_plot.png');
    fprintf('图像已保存到: results_nonmagnetic_plot.png\n');
end
