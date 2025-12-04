function [results] = calculate_NRW_with_gamma(s2p_file, d, options)
% 基于传播常数γ的NRW方法（使用公式 2-42 和 2-43）
%
% 输入参数:
%   s2p_file: S2P文件路径
%   d: 样品厚度 (m)
%   options: 可选参数
%       - lambda_c: 截止波长 (m)，默认 45.7e-3
%
% 核心公式:
%   (2-42) T = e^(-γd)
%   (2-43) γ = (2π/λ₀)√[εᵣμᵣ - (λ₀/λc)²]
%
% 反演流程:
%   S参数 → Γ → T → γ → (εᵣ, μᵣ)

    if nargin < 3
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;  % X波段
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    fprintf('========== NRW方法（基于γ公式）==========\n\n');
    fprintf('使用公式:\n');
    fprintf('  (2-42) T = exp(-γd)\n');
    fprintf('  (2-43) γ = (2π/λ₀)√[εᵣμᵣ - (λ₀/λc)²]\n\n');

    %% 读取S2P文件
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);

    lambda_0 = c ./ freq;

    %% 初始化
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    gamma_prop = zeros(n_points, 1);  % 传播常数
    mu_r = zeros(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);

    fprintf('开始计算...\n');

    for i = 1:n_points
        %% 步骤1: 从S参数计算反射系数Γ（公式 2-40, 2-41）
        X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));

        Gamma_p = X + sqrt(X^2 - 1);
        Gamma_m = X - sqrt(X^2 - 1);

        % 选择 |Γ| ≤ 1 的解
        if abs(Gamma_p) <= 1 && abs(Gamma_m) <= 1
            if abs(Gamma_m) < abs(Gamma_p)
                Gamma(i) = Gamma_m;
            else
                Gamma(i) = Gamma_p;
            end
        elseif abs(Gamma_p) <= 1
            Gamma(i) = Gamma_p;
        else
            Gamma(i) = Gamma_m;
        end

        %% 步骤2: 从S参数和Γ计算传输系数T
        % 标准NRW公式（这个不在你的图片里，但是必须的）
        T(i) = (S11(i) + S21(i) - Gamma(i)) / (1 - (S11(i) + S21(i)) * Gamma(i));

        %% 步骤3: 从T反推传播常数γ（公式 2-42 的逆运算）
        % T = e^(-γd)
        % γ = -ln(T) / d

        % 多分支相位解缠
        best_score = -inf;
        best_gamma = 0;

        for n_branch = -10:10
            % 考虑对数的多值性
            log_T = log(T(i)) + 2*pi*1j*n_branch;
            gamma_test = -log_T / d;

            %% 步骤4: 从γ计算εᵣμᵣ（公式 2-43 的逆运算）
            % γ = (2π/λ₀)√[εᵣμᵣ - (λ₀/λc)²]
            %
            % 设 K = γλ₀/(2π)
            % K² = εᵣμᵣ - (λ₀/λc)²
            % εᵣμᵣ = K² + (λ₀/λc)²

            K = gamma_test * lambda_0(i) / (2*pi);
            epsilon_mu_product = K^2 + (lambda_0(i) / lambda_c)^2;

            %% 步骤5: 从Γ和εᵣμᵣ分离出εᵣ和μᵣ
            % 使用阻抗关系（公式 2-44 到 2-47）
            %
            % Γ = (Z - Z₀) / (Z + Z₀)
            % Z/Z₀ = μᵣ/√(εᵣμᵣ - (λ₀/λc)²)
            % Z₀/Z₀ = 1/√(1 - (λ₀/λc)²)

            % 从Γ计算归一化阻抗
            % Γ = (z - 1)/(z + 1)  →  z = (1 + Γ)/(1 - Γ)
            z = (1 + Gamma(i)) / (1 - Gamma(i));

            % 根据传输线理论
            % z = √(μᵣ/εᵣ) / √(1 - (λ₀/λc)²)
            %
            % 所以 μᵣ/εᵣ = z² · (1 - (λ₀/λc)²)

            mu_epsilon_ratio = z^2 * (1 - (lambda_0(i) / lambda_c)^2);

            % 联立方程组：
            % εᵣμᵣ = epsilon_mu_product
            % μᵣ/εᵣ = mu_epsilon_ratio
            %
            % 解得：
            % μᵣ = √(epsilon_mu_product · mu_epsilon_ratio)
            % εᵣ = √(epsilon_mu_product / mu_epsilon_ratio)

            mu_r_test = sqrt(epsilon_mu_product * mu_epsilon_ratio);
            epsilon_r_test = sqrt(epsilon_mu_product / mu_epsilon_ratio);

            %% 评分（选择物理合理的分支）
            score = 0;

            % εᵣ实部必须为正
            if real(epsilon_r_test) > 0.5
                score = score + 20;
            end

            % μᵣ实部接近1（非磁性材料）
            if abs(real(mu_r_test) - 1) < 2
                score = score + 10;
            end

            % 虚部为正（损耗）
            if imag(epsilon_r_test) > 0 && imag(mu_r_test) > 0
                score = score + 10;
            end

            % 合理范围
            if real(epsilon_r_test) < 100 && abs(real(mu_r_test)) < 20
                score = score + 5;
            end

            % 连续性
            if i > 1 && ~isnan(epsilon_r(i-1))
                continuity = -abs(real(epsilon_r_test) - real(epsilon_r(i-1)));
                score = score + continuity * 0.1;
            end

            if score > best_score
                best_score = score;
                best_gamma = gamma_test;
                mu_r(i) = mu_r_test;
                epsilon_r(i) = epsilon_r_test;
            end
        end

        gamma_prop(i) = best_gamma;

        %% 计算损耗角正切
        if abs(real(epsilon_r(i))) > 1e-6
            tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
        else
            tan_delta(i) = NaN;
        end

        % 进度显示
        if mod(i, round(n_points/10)) == 0
            fprintf('  进度: %d/%d (%.1f%%)\n', i, n_points, i/n_points*100);
        end
    end

    fprintf('计算完成！\n\n');

    %% 验证：用计算出的εᵣ和μᵣ正向计算T，与测量值对比
    fprintf('========== 正向验证 ==========\n');
    T_reconstructed = zeros(n_points, 1);

    for i = 1:n_points
        % 用公式 (2-43) 计算γ
        term_inside = epsilon_r(i) * mu_r(i) - (lambda_0(i) / lambda_c)^2;
        gamma_calc = (2*pi / lambda_0(i)) * sqrt(term_inside);

        % 用公式 (2-42) 计算T
        T_reconstructed(i) = exp(-gamma_calc * d);
    end

    % 对比测量值和重建值
    T_error = abs(T - T_reconstructed);
    fprintf('T的重建误差:\n');
    fprintf('  平均误差: %.6f\n', mean(T_error));
    fprintf('  最大误差: %.6f\n', max(T_error));
    fprintf('  相对误差: %.2f%%\n\n', mean(T_error ./ abs(T)) * 100);

    %% 统计信息
    fprintf('========== 计算结果统计 ==========\n\n');

    valid_idx = real(epsilon_r) > 0 & real(epsilon_r) < 100;

    if sum(valid_idx) > 0
        fprintf('相对磁导率 μr (有效点 %d/%d):\n', sum(valid_idx), n_points);
        fprintf('  实部: %.4f ± %.4f\n', mean(real(mu_r(valid_idx))), std(real(mu_r(valid_idx))));

        fprintf('\n相对介电常数 εr (有效点 %d/%d):\n', sum(valid_idx), n_points);
        fprintf('  实部: %.4f ± %.4f\n', mean(real(epsilon_r(valid_idx))), std(real(epsilon_r(valid_idx))));

        fprintf('\n传播常数 γ:\n');
        fprintf('  衰减常数 α: %.4f Np/m\n', mean(real(gamma_prop(valid_idx))));
        fprintf('  相位常数 β: %.4f rad/m\n', mean(imag(gamma_prop(valid_idx))));
    end

    %% 保存结果
    output_file = 'results_NRW_gamma.xlsx';

    S11_dB = 20*log10(abs(S11));
    S21_dB = 20*log10(abs(S21));

    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(gamma_prop), imag(gamma_prop), abs(gamma_prop), ...
                         real(mu_r), imag(mu_r), abs(mu_r), ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, ...
                         abs(T_reconstructed), abs(T - T_reconstructed), ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'gamma_real', 'gamma_imag', 'gamma_mag', ...
                                          'mu_r_real', 'mu_r_imag', 'mu_r_mag', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_mag', ...
                                          'tan_delta', ...
                                          'T_reconstructed', 'T_error'});

    writetable(output_table, output_file);
    fprintf('\n结果已保存到: %s\n', output_file);

    %% 绘图
    plot_gamma_results(freq, S11_dB, S21_dB, Gamma, T, gamma_prop, mu_r, epsilon_r, tan_delta, T_reconstructed);

    %% 返回结果
    results.frequency = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.gamma = gamma_prop;
    results.mu_r = mu_r;
    results.epsilon_r = epsilon_r;
    results.tan_delta = tan_delta;
    results.T_reconstructed = T_reconstructed;

    fprintf('\n========== 完成 ==========\n');
end

function plot_gamma_results(freq, S11_dB, S21_dB, Gamma, T, gamma_prop, mu_r, epsilon_r, tan_delta, T_reconstructed)
    figure('Position', [100, 100, 1600, 1000], 'Name', 'NRW结果（基于γ公式）');

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

    % 子图2: Γ和T
    subplot(3, 3, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1.0, 'k--', 'LineWidth', 1);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|');
    title('反射/传输系数');

    % 子图3: T测量值vs重建值
    subplot(3, 3, 3);
    plot(freq/1e9, abs(T), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T_reconstructed), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|T|');
    legend('测量值', '重建值');
    title('T验证（公式2-42）');

    % 子图4: γ实部（衰减常数α）
    subplot(3, 3, 4);
    plot(freq/1e9, real(gamma_prop), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\alpha (Np/m)');
    title('衰减常数 α = Re(γ)');

    % 子图5: γ虚部（相位常数β）
    subplot(3, 3, 5);
    plot(freq/1e9, imag(gamma_prop), 'r-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\beta (rad/m)');
    title('相位常数 β = Im(γ)');

    % 子图6: |γ|
    subplot(3, 3, 6);
    plot(freq/1e9, abs(gamma_prop), 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('|\gamma| (1/m)');
    title('传播常数幅度');

    % 子图7: μr实部
    subplot(3, 3, 7);
    plot(freq/1e9, real(mu_r), 'b-', 'LineWidth', 2);
    hold on;
    yline(1.0, 'k--', '非磁性', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\mu_r''');
    title('相对磁导率（实部）');

    % 子图8: εr实部
    subplot(3, 3, 8);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 2);
    hold on;
    yline(0, 'r--', '物理下限', 'LineWidth', 0.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部）');

    % 子图9: tanδ
    subplot(3, 3, 9);
    valid_idx = abs(tan_delta) < 10;
    plot(freq(valid_idx)/1e9, tan_delta(valid_idx), 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');

    saveas(gcf, 'results_NRW_gamma_plot.png');
    fprintf('图像已保存到: results_NRW_gamma_plot.png\n');
end
