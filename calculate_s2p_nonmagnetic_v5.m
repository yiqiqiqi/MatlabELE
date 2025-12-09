function [results] = calculate_s2p_nonmagnetic_v5(s2p_file, d, options)
% 非磁性材料专用计算V5 - 添加自洽性验证
%
% 核心改进（V5）：
%   - 添加反向验证（Forward Validation）：从εr反推T，检查是否与测量匹配
%   - 自洽性好的k值获得巨大评分奖励
%   - 保持V4的公平评分和V3的相位连续性
%
% 参数：
%   s2p_file - S参数文件路径
%   d        - 样品厚度（米）
%   options  - 可选参数结构体
%       .lambda_c       - 截止波长（米），默认45.7mm
%       .k_search_range - k值搜索范围，默认20
%       .verbose        - 是否显示详细信息，默认true

    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;
    end

    if ~isfield(options, 'k_search_range')
        options.k_search_range = 20;
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    if options.verbose
        fprintf('========== 非磁性材料参数计算 V5 (自洽性验证) ==========\n\n');
        fprintf('核心改进: 反向验证机制（εr→T，检查与测量匹配度）\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  样品厚度 d = %.2f mm\n\n', d * 1000);
    end

    %% 读取S2P
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);
    lambda_0 = c ./ freq;

    %% 初始化输出
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    mu_r = ones(n_points, 1);
    tan_delta = zeros(n_points, 1);
    selected_k = zeros(n_points, 1);
    validation_error = zeros(n_points, 1);  % 记录验证误差

    %% 第一个点：全局搜索最佳k值（带自洽性检查）
    i = 1;
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

            % ===== 核心改进：反向验证 =====
            % 从εr反推γ和T，检查是否与测量的T匹配
            lambda_g_sq = lambda_0(i)^2 / epsilon_r_test - lambda_c^2;

            if real(lambda_g_sq) > 0
                % 传播模式
                gamma_reconstructed = 2*pi*1j / sqrt(lambda_g_sq);
            else
                % 截止模式（衰减）
                gamma_reconstructed = 2*pi / sqrt(-lambda_g_sq);
            end

            T_reconstructed = exp(-gamma_reconstructed * d);
            T_error = abs(T_reconstructed - T_test);

            % 评分（第一个点）
            score = 0;

            % 最高优先级：自洽性检查（Forward Validation）
            if T_error < 0.01
                score = score + 300;  % 极好的自洽性
            elseif T_error < 0.05
                score = score + 200;  % 好的自洽性
            elseif T_error < 0.1
                score = score + 100;  % 可接受的自洽性
            elseif T_error < 0.2
                score = score + 30;
            else
                score = score - T_error * 500;  % 不自洽，重罚
            end

            % εr'在物理合理范围
            eps_real = real(epsilon_r_test);

            if eps_real >= 0.8 && eps_real <= 20
                score = score + 100;
            elseif eps_real > 20 && eps_real <= 50
                score = score + 50;
            elseif eps_real >= 0.5 && eps_real < 0.8
                score = score + 30;
            else
                continue;  % 不合理，跳过
            end

            % εr''>=0 (损耗为正)
            eps_imag = imag(epsilon_r_test);
            if eps_imag >= -0.01
                score = score + 60;
                if eps_imag >= 0 && eps_imag <= 2
                    score = score + 20;
                end
            else
                score = score - abs(eps_imag) * 100;
            end

            % |T|<=1 (能量守恒)
            if abs(T_test) <= 1.0
                score = score + 80;
            elseif abs(T_test) <= 1.05
                score = score + 40;
            else
                score = score - (abs(T_test) - 1) * 200;
            end

            % |Γ|<=1 (能量守恒)
            if abs(Gamma_test) <= 1.0
                score = score + 60;
            elseif abs(Gamma_test) <= 1.05
                score = score + 20;
            else
                score = score - (abs(Gamma_test) - 1) * 150;
            end

            if score > best_score_global
                best_score_global = score;
                best_result_global.Gamma = Gamma_test;
                best_result_global.T = T_test;
                best_result_global.epsilon_r = epsilon_r_test;
                best_result_global.k = k;
                best_result_global.score = score;
                best_result_global.T_error = T_error;
            end
        end
    end

    % 应用第一个点的结果
    Gamma(1) = best_result_global.Gamma;
    T(1) = best_result_global.T;
    epsilon_r(1) = best_result_global.epsilon_r;
    selected_k(1) = best_result_global.k;
    validation_error(1) = best_result_global.T_error;

    if options.verbose
        fprintf('  第一个点 (%.2f GHz):\n', freq(1)/1e9);
        fprintf('    选择: k=%d, εr''=%.3f\n', selected_k(1), real(epsilon_r(1)));
        fprintf('    验证误差: |T_rec - T_mea| = %.4f\n', validation_error(1));
        fprintf('    得分: %.0f\n', best_result_global.score);
    end

    % 计算tanδ
    if abs(real(epsilon_r(1))) > 1e-6
        tan_delta(1) = imag(epsilon_r(1)) / real(epsilon_r(1));
    end

    %% 后续点：优先使用前一个点的k值（相位连续性）+ 自洽性检查
    for i = 2:n_points
        try
            X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
            Gamma_p = X + sqrt(X^2 - 1);
            Gamma_m = X - sqrt(X^2 - 1);

            gamma_candidates = [Gamma_p; Gamma_m];

            prev_k = selected_k(i-1);
            k_priority = [prev_k, prev_k-1, prev_k+1];

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

                    % ===== 反向验证 =====
                    lambda_g_sq = lambda_0(i)^2 / epsilon_r_test - lambda_c^2;

                    if real(lambda_g_sq) > 0
                        gamma_reconstructed = 2*pi*1j / sqrt(lambda_g_sq);
                    else
                        gamma_reconstructed = 2*pi / sqrt(-lambda_g_sq);
                    end

                    T_reconstructed = exp(-gamma_reconstructed * d);
                    T_error = abs(T_reconstructed - T_test);

                    %% 评分（强调连续性+自洽性）
                    score = 0;

                    % 最高优先级1：自洽性
                    if T_error < 0.01
                        score = score + 250;
                    elseif T_error < 0.05
                        score = score + 180;
                    elseif T_error < 0.1
                        score = score + 100;
                    elseif T_error < 0.2
                        score = score + 30;
                    else
                        score = score - T_error * 400;
                    end

                    % 最高优先级2：k值连续性
                    if k == prev_k
                        score = score + 250;
                    elseif abs(k - prev_k) == 1
                        score = score + 80;
                    else
                        score = score - abs(k - prev_k) * 150;
                    end

                    % εr'连续性
                    eps_real = real(epsilon_r_test);
                    prev_eps_real = real(epsilon_r(i-1));
                    diff_eps = abs(eps_real - prev_eps_real);

                    if diff_eps < 0.5
                        score = score + 120;
                    elseif diff_eps < 1
                        score = score + 80;
                    elseif diff_eps < 2
                        score = score + 40;
                    elseif diff_eps < 5
                        score = score + 10;
                    else
                        score = score - diff_eps * 30;
                    end

                    % εr'在合理范围
                    if eps_real >= 0.8 && eps_real <= 100
                        score = score + 40;
                    else
                        score = score - 80;
                    end

                    % εr''>=0
                    eps_imag = imag(epsilon_r_test);
                    if eps_imag >= -0.01
                        score = score + 30;
                    else
                        score = score - abs(eps_imag) * 80;
                    end

                    % |T|<=1
                    if abs(T_test) <= 1.0
                        score = score + 50;
                    elseif abs(T_test) <= 1.05
                        score = score + 20;
                    else
                        score = score - (abs(T_test) - 1) * 150;
                    end

                    % |Γ|<=1
                    if abs(Gamma_test) <= 1.0
                        score = score + 30;
                    elseif abs(Gamma_test) <= 1.05
                        score = score + 10;
                    end

                    if score > best_score
                        best_score = score;
                        best_result.Gamma = Gamma_test;
                        best_result.T = T_test;
                        best_result.epsilon_r = epsilon_r_test;
                        best_result.k = k;
                        best_result.T_error = T_error;
                    end
                end
            end

            %% 策略2: 如果优先k值都不好，扩大搜索
            if best_score < -50
                for g_idx = 1:length(gamma_candidates)
                    Gamma_test = gamma_candidates(g_idx);
                    T_test = (S11(i) + S21(i) - Gamma_test) / ...
                             (1 - (S11(i) + S21(i)) * Gamma_test);

                    for k = (prev_k-5):(prev_k+5)
                        if ismember(k, k_priority)
                            continue;
                        end

                        lgT = log(1/T_test) + 2*pi*1j*k;
                        inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
                        epsilon_r_test = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

                        % 反向验证
                        lambda_g_sq = lambda_0(i)^2 / epsilon_r_test - lambda_c^2;

                        if real(lambda_g_sq) > 0
                            gamma_reconstructed = 2*pi*1j / sqrt(lambda_g_sq);
                        else
                            gamma_reconstructed = 2*pi / sqrt(-lambda_g_sq);
                        end

                        T_reconstructed = exp(-gamma_reconstructed * d);
                        T_error = abs(T_reconstructed - T_test);

                        score = 0;

                        % 自洽性
                        if T_error < 0.05
                            score = score + 150;
                        elseif T_error < 0.1
                            score = score + 80;
                        else
                            score = score - T_error * 300;
                        end

                        % k跳变惩罚
                        score = score - abs(k - prev_k) * 150;

                        % εr'连续性
                        diff_eps = abs(real(epsilon_r_test) - real(epsilon_r(i-1)));
                        if diff_eps < 1
                            score = score + 80;
                        elseif diff_eps < 3
                            score = score + 40;
                        else
                            score = score - diff_eps * 20;
                        end

                        % 物理约束
                        if real(epsilon_r_test) >= 0.8 && real(epsilon_r_test) <= 100
                            score = score + 40;
                        end

                        if imag(epsilon_r_test) >= -0.01
                            score = score + 20;
                        end

                        if abs(T_test) <= 1.05
                            score = score + 30;
                        end

                        if abs(Gamma_test) <= 1.05
                            score = score + 15;
                        end

                        if score > best_score
                            best_score = score;
                            best_result.Gamma = Gamma_test;
                            best_result.T = T_test;
                            best_result.epsilon_r = epsilon_r_test;
                            best_result.k = k;
                            best_result.T_error = T_error;
                        end
                    end
                end
            end

            % 应用结果
            if ~isempty(fieldnames(best_result))
                Gamma(i) = best_result.Gamma;
                T(i) = best_result.T;
                epsilon_r(i) = best_result.epsilon_r;
                selected_k(i) = best_result.k;
                validation_error(i) = best_result.T_error;
            else
                % 极端情况
                Gamma(i) = Gamma(i-1);
                T(i) = T(i-1);
                epsilon_r(i) = epsilon_r(i-1);
                selected_k(i) = selected_k(i-1);
                validation_error(i) = validation_error(i-1);
            end

            % 计算tanδ
            if abs(real(epsilon_r(i))) > 1e-6
                tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
            end

        catch ME
            warning('频点 %d 计算失败: %s', i, ME.message);
            if i > 1
                Gamma(i) = Gamma(i-1);
                T(i) = T(i-1);
                epsilon_r(i) = epsilon_r(i-1);
                selected_k(i) = selected_k(i-1);
                tan_delta(i) = tan_delta(i-1);
                validation_error(i) = validation_error(i-1);
            end
        end
    end

    %% 诊断输出
    if options.verbose
        fprintf('\n========== 诊断信息 ==========\n');

        % k值跳变
        k_jumps = diff(selected_k);
        n_jumps = sum(k_jumps ~= 0);
        fprintf('  k值跳变次数: %d / %d 个频点\n', n_jumps, n_points-1);

        if n_jumps > 0
            jump_indices = find(k_jumps ~= 0);
            fprintf('  跳变位置:\n');
            for j = 1:min(5, length(jump_indices))
                idx = jump_indices(j);
                fprintf('    频点 %d (%.2f GHz): k=%d -> k=%d, Δεr''=%.3f\n', ...
                        idx, freq(idx)/1e9, selected_k(idx), selected_k(idx+1), ...
                        real(epsilon_r(idx+1)) - real(epsilon_r(idx)));
            end
            if length(jump_indices) > 5
                fprintf('    ... 还有 %d 个跳变点\n', length(jump_indices) - 5);
            end
        end

        % εr'统计
        fprintf('\n  εr''范围: %.3f ~ %.3f\n', min(real(epsilon_r)), max(real(epsilon_r)));
        fprintf('  εr''均值: %.3f\n', mean(real(epsilon_r)));
        fprintf('  εr''标准差: %.4f\n', std(real(epsilon_r)));

        % 验证误差
        fprintf('\n  验证误差统计:\n');
        fprintf('    平均误差: %.4f\n', mean(validation_error));
        fprintf('    最大误差: %.4f\n', max(validation_error));
        fprintf('    误差<0.01的点: %d / %d (%.1f%%)\n', ...
                sum(validation_error < 0.01), n_points, ...
                100*sum(validation_error < 0.01)/n_points);
        fprintf('========================================\n');
    end

    %% 封装结果
    results = struct();
    results.freq = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.epsilon_r = epsilon_r;
    results.mu_r = mu_r;
    results.tan_delta = tan_delta;
    results.lambda_c = lambda_c;
    results.d = d;
    results.selected_k = selected_k;
    results.validation_error = validation_error;

    %% 绘图
    figure('Name', 'S参数与材料参数 (V5自洽性验证)', 'Position', [100, 100, 1400, 900]);

    % S参数 (dB)
    subplot(3, 3, 1);
    plot(freq/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, 20*log10(abs(S21)), 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    title('S参数');
    legend('S_{11}', 'S_{21}');

    % Γ和T
    subplot(3, 3, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1, 'k--', '理想上限');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('模值');
    title('反射/传输系数');
    legend('|\Gamma|', '|T|');

    % 相位分支k值
    subplot(3, 3, 3);
    plot(freq/1e9, selected_k, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('k值 (相位分支)');
    title('相位分支演化');

    % εr' (实部)
    subplot(3, 3, 4);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 2);
    hold on;
    if min(real(epsilon_r)) < 2
        yline(1, 'k--', '空气');
    end
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon''_r');
    title('相对介电常数（实部）');

    % εr'' (虚部)
    subplot(3, 3, 5);
    plot(freq/1e9, imag(epsilon_r), 'r-', 'LineWidth', 1.5);
    hold on;
    yline(0, 'k--');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon''''_r');
    title('相对介电常数（虚部）');

    % tanδ
    subplot(3, 3, 6);
    plot(freq/1e9, tan_delta, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    % 验证误差（关键！）
    subplot(3, 3, 7);
    semilogy(freq/1e9, validation_error, 'r-', 'LineWidth', 1.5);
    hold on;
    yline(0.01, 'g--', '良好阈值');
    yline(0.05, 'y--', '可接受阈值');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|T_{重构} - T_{测量}|');
    title('自洽性验证误差（越小越好）');

    % k值与εr'对比
    subplot(3, 3, 8);
    yyaxis left;
    plot(freq/1e9, selected_k, 'b-', 'LineWidth', 2);
    ylabel('k值', 'Color', 'b');

    yyaxis right;
    plot(freq/1e9, real(epsilon_r), 'r-', 'LineWidth', 1.5);
    ylabel('\epsilon''_r', 'Color', 'r');

    grid on;
    xlabel('频率 (GHz)');
    title('k值与εr''对比');

    % 验证误差 vs k值（散点图）
    subplot(3, 3, 9);
    scatter(selected_k, validation_error, 50, freq/1e9, 'filled');
    xlabel('k值');
    ylabel('验证误差');
    title('验证误差 vs k值（颜色=频率GHz）');
    colorbar;
    grid on;
    set(gca, 'YScale', 'log');

end
