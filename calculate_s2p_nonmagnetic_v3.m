function [results] = calculate_s2p_nonmagnetic_v3(s2p_file, d, options)
% 非磁性材料专用计算V3 - 强制相位连续性
%
% 核心改进：
%   - 强制k值连续性（防止相位分支跳变导致εr突变）
%   - 初始化策略：第一个点careful选择，后续沿着同一分支
%   - 只在必要时允许k值±1跳变

    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;
    end

    if ~isfield(options, 'k_search_range')
        options.k_search_range = 20;  % 初始搜索范围
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    if options.verbose
        fprintf('========== 非磁性材料参数计算 V3 ==========\n\n');
        fprintf('核心改进: 强制相位连续性，防止k值跳变\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  样品厚度 d = %.2f mm\n\n', d * 1000);
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

    fprintf('开始计算...\n');

    %% 第一个点：全局搜索找最优k值
    i = 1;
    fprintf('  初始化：全局搜索第一个点的最佳k值...\n');

    X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
    Gamma_p = X + sqrt(X^2 - 1);
    Gamma_m = X - sqrt(X^2 - 1);

    gamma_candidates = [Gamma_p; Gamma_m];

    best_score_global = -inf;
    best_result_global = struct();

    for g_idx = 1:length(gamma_candidates)
        Gamma_test = gamma_candidates(g_idx);
        T_test = (S11(i) + S21(i) - Gamma_test) / ...
                 (1 - (S11(i) + S21(i)) * Gamma_test);

        for k = -options.k_search_range:options.k_search_range
            lgT = log(1/T_test) + 2*pi*1j*k;
            inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
            epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

            % 评分（第一个点）
            score = 0;

            % εr'在合理范围（2-50）
            eps_real = real(epsilon_r_test);
            if eps_real >= 2 && eps_real <= 50
                score = score + 100;
            elseif eps_real >= 1 && eps_real < 2
                score = score + 50;
            elseif eps_real > 50 && eps_real <= 100
                score = score + 30;
            else
                continue;  % 不合理，跳过
            end

            % εr''>=0
            eps_imag = imag(epsilon_r_test);
            if eps_imag >= 0
                score = score + 50;
            else
                score = score - abs(eps_imag) * 100;
            end

            % |T|<=1
            if abs(T_test) <= 1.0
                score = score + 50;
            elseif abs(T_test) <= 1.1
                score = score + 20;
            else
                score = score - (abs(T_test) - 1) * 100;
            end

            % |Γ|<=1
            if abs(Gamma_test) <= 1.0
                score = score + 30;
            elseif abs(Gamma_test) <= 1.1
                score = score + 10;
            end

            if score > best_score_global
                best_score_global = score;
                best_result_global.Gamma = Gamma_test;
                best_result_global.T = T_test;
                best_result_global.epsilon_r = epsilon_r_test;
                best_result_global.k = k;
            end
        end
    end

    % 应用第一个点的结果
    Gamma(1) = best_result_global.Gamma;
    T(1) = best_result_global.T;
    epsilon_r(1) = best_result_global.epsilon_r;
    selected_k(1) = best_result_global.k;

    fprintf('  第一个点选择: k=%d, εr''=%.2f\n', selected_k(1), real(epsilon_r(1)));

    % 计算tanδ
    if abs(real(epsilon_r(1))) > 1e-6
        tan_delta(1) = imag(epsilon_r(1)) / real(epsilon_r(1));
    end

    %% 后续点：优先使用前一个点的k值（相位连续性）
    for i = 2:n_points
        try
            X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
            Gamma_p = X + sqrt(X^2 - 1);
            Gamma_m = X - sqrt(X^2 - 1);

            gamma_candidates = [Gamma_p; Gamma_m];

            % 当前点优先搜索范围：前一个k值 ± 1
            prev_k = selected_k(i-1);
            k_priority = [prev_k, prev_k-1, prev_k+1];  % 优先级顺序

            best_score = -inf;
            best_result = struct();

            %% 策略1: 优先尝试相同k值和相邻k值
            for g_idx = 1:length(gamma_candidates)
                Gamma_test = gamma_candidates(g_idx);
                T_test = (S11(i) + S21(i) - Gamma_test) / ...
                         (1 - (S11(i) + S21(i)) * Gamma_test);

                for k_idx = 1:length(k_priority)
                    k = k_priority(k_idx);

                    lgT = log(1/T_test) + 2*pi*1j*k;
                    inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
                    epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

                    %% 评分（强调连续性）
                    score = 0;

                    % 最高优先级：k值连续性
                    if k == prev_k
                        score = score + 200;  // 巨大奖励！
                    elseif abs(k - prev_k) == 1
                        score = score + 50;   // 允许±1跳变
                    else
                        score = score - abs(k - prev_k) * 100;  // 惩罚大跳变
                    end

                    % εr'连续性
                    eps_real = real(epsilon_r_test);
                    prev_eps_real = real(epsilon_r(i-1));
                    diff_eps = abs(eps_real - prev_eps_real);

                    if diff_eps < 1
                        score = score + 100;
                    elseif diff_eps < 3
                        score = score + 50;
                    elseif diff_eps < 5
                        score = score + 20;
                    else
                        score = score - diff_eps * 20;
                    end

                    % εr'合理范围
                    if eps_real >= 1 && eps_real <= 100
                        score = score + 50;
                    else
                        score = score - 100;
                    end

                    % εr''>=0
                    eps_imag = imag(epsilon_r_test);
                    if eps_imag >= 0
                        score = score + 30;
                    else
                        score = score - abs(eps_imag) * 50;
                    end

                    % |T|<=1
                    if abs(T_test) <= 1.0
                        score = score + 40;
                    elseif abs(T_test) <= 1.1
                        score = score + 15;
                    else
                        score = score - (abs(T_test) - 1) * 80;
                    end

                    % |Γ|<=1
                    if abs(Gamma_test) <= 1.0
                        score = score + 20;
                    elseif abs(Gamma_test) <= 1.1
                        score = score + 5;
                    end

                    if score > best_score
                        best_score = score;
                        best_result.Gamma = Gamma_test;
                        best_result.T = T_test;
                        best_result.epsilon_r = epsilon_r_test;
                        best_result.k = k;
                    end
                end
            end

            %% 策略2: 如果优先k值都不好，扩大搜索（但仍然惩罚跳变）
            if best_score < -50  % 阈值：如果分数太低
                for g_idx = 1:length(gamma_candidates)
                    Gamma_test = gamma_candidates(g_idx);
                    T_test = (S11(i) + S21(i) - Gamma_test) / ...
                             (1 - (S11(i) + S21(i)) * Gamma_test);

                    for k = (prev_k-5):(prev_k+5)  % 扩大到±5
                        if ismember(k, k_priority)
                            continue;  % 已经测试过了
                        end

                        lgT = log(1/T_test) + 2*pi*1j*k;
                        inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
                        epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

                        % 同样的评分，但k值跳变惩罚更大
                        score = -abs(k - prev_k) * 150;  // 强烈惩罚

                        eps_real = real(epsilon_r_test);
                        if eps_real >= 1 && eps_real <= 100
                            score = score + 50;
                        else
                            continue;
                        end

                        if imag(epsilon_r_test) >= 0
                            score = score + 30;
                        end

                        if abs(T_test) <= 1.1
                            score = score + 30;
                        end

                        if score > best_score
                            best_score = score;
                            best_result.Gamma = Gamma_test;
                            best_result.T = T_test;
                            best_result.epsilon_r = epsilon_r_test;
                            best_result.k = k;
                        end
                    end
                end
            end

            %% 应用最佳结果
            if ~isempty(fieldnames(best_result))
                Gamma(i) = best_result.Gamma;
                T(i) = best_result.T;
                epsilon_r(i) = best_result.epsilon_r;
                selected_k(i) = best_result.k;

                % 记录k值跳变
                if selected_k(i) ~= prev_k
                    warnings_list{i} = sprintf('k跳变: %d→%d', prev_k, selected_k(i));
                end

                % 其他警告
                if abs(T(i)) > 1.05
                    if isempty(warnings_list{i})
                        warnings_list{i} = sprintf('|T|=%.3f', abs(T(i)));
                    else
                        warnings_list{i} = sprintf('%s; |T|=%.3f', warnings_list{i}, abs(T(i)));
                    end
                end
            else
                warnings_list{i} = '无合理解';
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
            fprintf('  进度: %d/%d (%.1f%%), 当前k=%d\n', i, n_points, i/n_points*100, selected_k(i));
        end
    end

    fprintf('计算完成！\n\n');

    %% 诊断报告
    if options.verbose
        fprintf('========== 诊断报告 ==========\n\n');

        % k值跳变统计
        k_jumps = abs(diff(selected_k));
        n_jumps = sum(k_jumps > 0);
        n_big_jumps = sum(k_jumps > 1);

        fprintf('相位连续性分析:\n');
        fprintf('  k值范围: %d 到 %d\n', min(selected_k), max(selected_k));
        fprintf('  k值跳变次数: %d (%.1f%%)\n', n_jumps, n_jumps/(n_points-1)*100);
        fprintf('  大跳变(>1): %d (%.1f%%)\n', n_big_jumps, n_big_jumps/(n_points-1)*100);
        fprintf('  平均k值: %.2f\n\n', mean(selected_k));

        % 物理约束
        valid_T = sum(abs(T) <= 1.05);
        valid_eps = sum(real(epsilon_r) >= 1 & real(epsilon_r) < 100);
        valid_eps_imag = sum(imag(epsilon_r) >= 0);

        fprintf('物理约束符合率:\n');
        fprintf('  |T| ≤ 1.05: %d/%d (%.1f%%)\n', valid_T, n_points, valid_T/n_points*100);
        fprintf('  1 ≤ εr'' < 100: %d/%d (%.1f%%)\n', valid_eps, n_points, valid_eps/n_points*100);
        fprintf('  εr'''' ≥ 0: %d/%d (%.1f%%)\n\n', valid_eps_imag, n_points, valid_eps_imag/n_points*100);

        % 材料参数统计
        good_idx = (abs(T) <= 1.05) & (real(epsilon_r) >= 1) & (imag(epsilon_r) >= 0);
        if sum(good_idx) > 0
            fprintf('高质量数据点: %d/%d (%.1f%%)\n', sum(good_idx), n_points, sum(good_idx)/n_points*100);
            fprintf('推荐频率范围: %.3f - %.3f GHz\n\n', ...
                    min(freq(good_idx))/1e9, max(freq(good_idx))/1e9);

            fprintf('材料参数（高质量点）:\n');
            fprintf('  εr'' = %.3f ± %.3f\n', mean(real(epsilon_r(good_idx))), std(real(epsilon_r(good_idx))));
            fprintf('  εr'''' = %.3f ± %.3f\n', mean(imag(epsilon_r(good_idx))), std(imag(epsilon_r(good_idx))));
            fprintf('  tanδ = %.4f ± %.4f\n', mean(tan_delta(good_idx)), std(tan_delta(good_idx)));
        end
    end

    %% 保存结果
    S11_dB = 20*log10(abs(S11));
    S21_dB = 20*log10(abs(S21));

    output_file = 'results_nonmagnetic_v3.xlsx';
    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, selected_k, warnings_list, ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_mag', ...
                                          'tan_delta', 'Selected_k', 'Warnings'});

    writetable(output_table, output_file);
    fprintf('\n结果已保存到: %s\n', output_file);

    %% 绘图
    plot_results_v3(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, selected_k, good_idx);

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
    results.warnings = warnings_list;
    results.good_idx = good_idx;

    fprintf('\n========== 完成 ==========\n');
end

function plot_results_v3(freq, S11_dB, S21_dB, Gamma, T, epsilon_r, tan_delta, selected_k, good_idx)
    figure('Position', [100, 100, 1600, 1000], 'Name', '非磁性材料分析 V3 - 相位连续性');

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
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|', 'Location', 'best');
    title('反射/传输系数');

    % 子图3: k值演化（关键！）
    subplot(3, 3, 3);
    plot(freq/1e9, selected_k, 'k-', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('k值（相位分支）');
    title('相位分支演化（应该连续）');

    % 标记跳变点
    k_jumps_idx = find(abs(diff(selected_k)) > 0);
    if ~isempty(k_jumps_idx)
        hold on;
        plot(freq(k_jumps_idx+1)/1e9, selected_k(k_jumps_idx+1), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        legend('k值', 'k跳变点', 'Location', 'best');
    end

    % 子图4: εr实部
    subplot(3, 3, 4);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 2);
    hold on;
    yline(1.0, 'r--', '物理下限', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, real(epsilon_r(good_idx)), 'g.', 'MarkerSize', 8);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部）- 应平滑');

    % 子图5: εr虚部
    subplot(3, 3, 5);
    plot(freq/1e9, imag(epsilon_r), 'r-', 'LineWidth', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 0.5);
    if exist('good_idx', 'var')
        plot(freq(good_idx)/1e9, imag(epsilon_r(good_idx)), 'g.', 'MarkerSize', 8);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''''');
    title('相对介电常数（虚部）');

    % 子图6: tanδ
    subplot(3, 3, 6);
    valid_tan = abs(tan_delta) < 5;
    plot(freq(valid_tan)/1e9, tan_delta(valid_tan), 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    % 子图7: k值与εr'的关系
    subplot(3, 3, 7);
    yyaxis left
    plot(freq/1e9, selected_k, 'k-', 'LineWidth', 1.5);
    ylabel('k值');
    yyaxis right
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 1.5);
    ylabel('\epsilon_r''');
    grid on;
    xlabel('频率 (GHz)');
    title('k值与εr''对比（查看相关性）');

    % 子图8: εr'的变化率
    subplot(3, 3, 8);
    deps_df = [0; diff(real(epsilon_r))];
    plot(freq/1e9, abs(deps_df), 'b-', 'LineWidth', 1.5);
    hold on;
    % 标记大跳变点
    big_change = abs(deps_df) > 5;
    if sum(big_change) > 0
        plot(freq(big_change)/1e9, abs(deps_df(big_change)), 'ro', 'MarkerSize', 8);
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|Δεr''|');
    title('εr''变化率（大跳变=红圈）');
    legend('变化率', '大跳变', 'Location', 'best');

    % 子图9: k跳变与εr'跳变的关联
    subplot(3, 3, 9);
    k_change = [0; abs(diff(selected_k))];
    scatter(k_change, abs(deps_df), 50, freq/1e9, 'filled');
    colorbar;
    xlabel('k值变化');
    ylabel('|Δεr''|');
    title('k跳变 vs εr''跳变（颜色=频率GHz）');
    grid on;

    saveas(gcf, 'results_nonmagnetic_v3_plot.png');
    fprintf('图像已保存到: results_nonmagnetic_v3_plot.png\n');
end
