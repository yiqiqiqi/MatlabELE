function [results] = calculate_s2p_diagnostic(s2p_file, d, k_values, options)
% 诊断版本：显示不同k值对应的εr'计算结果，不做自动选择
%
% 目的：让用户直观看到不同相位分支k的影响，自行判断哪个k值合理
%
% 参数：
%   s2p_file - S参数文件路径
%   d        - 样品厚度（米）
%   k_values - 要尝试的k值数组，默认 [-2, -1, 0, 1, 2]
%   options  - 可选参数
%       .lambda_c - 截止波长（米），默认45.7mm
%       .verbose  - 是否显示详细信息，默认true

    if nargin < 3 || isempty(k_values)
        k_values = [-2, -1, 0, 1, 2];
    end

    if nargin < 4
        options = struct();
    end

    if ~isfield(options, 'lambda_c')
        options.lambda_c = 45.7e-3;
    end

    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    %% 物理常数
    c = 299792458;
    lambda_c = options.lambda_c;

    if options.verbose
        fprintf('========== S参数诊断分析（无自动选择） ==========\n\n');
        fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
        fprintf('  样品厚度 d = %.2f mm\n', d * 1000);
        fprintf('  尝试的k值: [%s]\n', num2str(k_values));
        fprintf('  目的：展示不同k值的εr''，由用户判断哪个合理\n\n');
    end

    %% 读取S2P
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);
    n_points = length(freq);
    lambda_0 = c ./ freq;

    %% 初始化
    n_k = length(k_values);

    % 对每个k值，存储εr'
    epsilon_r_real_all = zeros(n_points, n_k);

    % 共用的Γ和T
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);

    %% 对每个频率点计算
    for i = 1:n_points
        % 计算Γ（选择|Γ|≤1的）
        X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));
        Gamma_p = X + sqrt(X^2 - 1);
        Gamma_m = X - sqrt(X^2 - 1);

        % 选择满足|Γ|≤1的，优先选模较小的
        if abs(Gamma_p) <= 1.0 && abs(Gamma_m) <= 1.0
            if abs(Gamma_p) < abs(Gamma_m)
                Gamma(i) = Gamma_p;
            else
                Gamma(i) = Gamma_m;
            end
        elseif abs(Gamma_p) <= 1.0
            Gamma(i) = Gamma_p;
        elseif abs(Gamma_m) <= 1.0
            Gamma(i) = Gamma_m;
        else
            % 都不满足，选模较小的
            if abs(Gamma_p) < abs(Gamma_m)
                Gamma(i) = Gamma_p;
            else
                Gamma(i) = Gamma_m;
            end
        end

        % 计算T
        T(i) = (S11(i) + S21(i) - Gamma(i)) / ...
               (1 - (S11(i) + S21(i)) * Gamma(i));

        % 对每个k值计算εr'
        for k_idx = 1:n_k
            k = k_values(k_idx);

            lgT = log(1/T(i)) + 2*pi*1j*k;
            inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
            epsilon_r = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

            epsilon_r_real_all(i, k_idx) = real(epsilon_r);
        end
    end

    %% 统计信息
    if options.verbose
        fprintf('========== 不同k值的εr''统计 ==========\n');
        for k_idx = 1:n_k
            k = k_values(k_idx);
            eps_r = epsilon_r_real_all(:, k_idx);
            fprintf('  k=%+2d: εr''范围 [%.2f, %.2f], 均值=%.2f, 标准差=%.3f\n', ...
                    k, min(eps_r), max(eps_r), mean(eps_r), std(eps_r));
        end
        fprintf('\n提示：\n');
        fprintf('  - 空气应该εr''≈1\n');
        fprintf('  - 介电材料εr''>1且应该平滑变化\n');
        fprintf('  - 选择最符合物理预期的k值\n');
        fprintf('==========================================\n');
    end

    %% 封装结果
    results = struct();
    results.freq = freq;
    results.S11 = S11;
    results.S21 = S21;
    results.Gamma = Gamma;
    results.T = T;
    results.k_values = k_values;
    results.epsilon_r_real_all = epsilon_r_real_all;
    results.lambda_c = lambda_c;
    results.d = d;

    %% 绘图
    figure('Name', 'S参数诊断分析（不同k值对比）', 'Position', [100, 100, 1200, 800]);

    % S参数 (dB)
    subplot(2, 3, 1);
    plot(freq/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, 20*log10(abs(S21)), 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    title('S参数');
    legend('S_{11}', 'S_{21}');

    % Γ和T的模
    subplot(2, 3, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5);
    yline(1, 'k--', '理想上限');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('模值');
    title('反射/传输系数');
    legend('|\Gamma|', '|T|', 'Location', 'best');

    % Γ的相位
    subplot(2, 3, 3);
    plot(freq/1e9, angle(Gamma)*180/pi, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相位 (度)');
    title('\Gamma的相位');

    % T的相位
    subplot(2, 3, 4);
    plot(freq/1e9, angle(T)*180/pi, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相位 (度)');
    title('T的相位');

    % 核心：不同k值的εr'（大图）
    subplot(2, 3, [5, 6]);
    colors = lines(n_k);
    hold on;

    for k_idx = 1:n_k
        k = k_values(k_idx);
        plot(freq/1e9, epsilon_r_real_all(:, k_idx), ...
             'LineWidth', 2, 'Color', colors(k_idx, :), ...
             'DisplayName', sprintf('k=%+d', k));
    end

    % 参考线
    yline(1, 'k--', '空气', 'LineWidth', 1, 'Alpha', 0.5);

    grid on;
    xlabel('频率 (GHz)', 'FontSize', 12);
    ylabel('\epsilon''_r', 'FontSize', 12);
    title('不同k值对应的相对介电常数实部（选择平滑且符合预期的曲线）', 'FontSize', 13);
    legend('Location', 'best', 'FontSize', 11);
    ylim([min(epsilon_r_real_all(:))-1, max(epsilon_r_real_all(:))+1]);

    % 在命令行输出提示
    fprintf('\n请查看图中"不同k值对应的相对介电常数实部"：\n');
    fprintf('  - 选择物理合理的曲线（平滑、符合材料特性）\n');
    fprintf('  - 空气：应选择εr''≈1的k值\n');
    fprintf('  - 介电材料：应选择给出稳定εr''>1的k值\n\n');

end
