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

    % 对每个k值，存储εr'和μr'
    epsilon_r_real_all = zeros(n_points, n_k);
    mu_r_real_all = zeros(n_points, n_k);

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

        % 计算阻抗比（从Γ）
        z = (1 + Gamma(i)) / (1 - Gamma(i));

        % 对每个k值计算εr和μr
        for k_idx = 1:n_k
            k = k_values(k_idx);

            % 计算εr（从T和k）
            lgT = log(1/T(i)) + 2*pi*1j*k;
            inv_lambda_sq = -((1/(2*pi*d)) * lgT)^2;
            epsilon_r = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

            % 计算μr（从εr和z）
            % 波导中: z = sqrt(μr/εr) / sqrt(1-(λ0/λc)²)
            % 所以: μr/εr = z² * (1-(λ0/λc)²)
            lambda_ratio_sq = (lambda_0(i)/lambda_c)^2;
            if lambda_ratio_sq < 1
                mu_epsilon_ratio = z^2 * (1 - lambda_ratio_sq);
                mu_r = sqrt(epsilon_r * mu_epsilon_ratio);
            else
                % 截止模式
                mu_r = 1 + 0i;  % 默认值
            end

            epsilon_r_real_all(i, k_idx) = real(epsilon_r);
            mu_r_real_all(i, k_idx) = real(mu_r);
        end
    end

    %% 统计信息
    if options.verbose
        fprintf('========== 不同k值的εr''和μr''统计 ==========\n');
        for k_idx = 1:n_k
            k = k_values(k_idx);
            eps_r = epsilon_r_real_all(:, k_idx);
            mu_r = mu_r_real_all(:, k_idx);
            fprintf('  k=%+2d: εr''∈[%.2f, %.2f] (均值%.2f), μr''∈[%.2f, %.2f] (均值%.2f)\n', ...
                    k, min(eps_r), max(eps_r), mean(eps_r), ...
                    min(mu_r), max(mu_r), mean(mu_r));
        end
        fprintf('\n提示：\n');
        fprintf('  - 空气：εr''≈1, μr''≈1\n');
        fprintf('  - 非磁性介电材料：εr''>1, μr''≈1\n');
        fprintf('  - 选择给出物理合理结果的k值\n');
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
    results.mu_r_real_all = mu_r_real_all;
    results.lambda_c = lambda_c;
    results.d = d;

    %% 绘图
    figure('Name', 'S参数诊断分析（不同k值对比）', 'Position', [100, 100, 1400, 900]);

    % S参数 (dB)
    subplot(3, 2, 1);
    plot(freq/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(freq/1e9, 20*log10(abs(S21)), 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    title('S参数');
    legend('S_{11}', 'S_{21}');

    % Γ和T的模
    subplot(3, 2, 2);
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
    subplot(3, 2, 3);
    plot(freq/1e9, angle(Gamma)*180/pi, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相位 (度)');
    title('\Gamma的相位');

    % T的相位
    subplot(3, 2, 4);
    plot(freq/1e9, angle(T)*180/pi, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相位 (度)');
    title('T的相位');

    % 核心1：不同k值的εr'
    subplot(3, 2, 5);
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
    xlabel('频率 (GHz)', 'FontSize', 11);
    ylabel('\epsilon''_r', 'FontSize', 11);
    title('不同k值对应的相对介电常数实部', 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 10);
    ylim([min(epsilon_r_real_all(:))-1, max(epsilon_r_real_all(:))+1]);

    % 核心2：不同k值的μr'
    subplot(3, 2, 6);
    hold on;

    for k_idx = 1:n_k
        k = k_values(k_idx);
        plot(freq/1e9, mu_r_real_all(:, k_idx), ...
             'LineWidth', 2, 'Color', colors(k_idx, :), ...
             'DisplayName', sprintf('k=%+d', k));
    end

    % 参考线
    yline(1, 'k--', '非磁性材料', 'LineWidth', 1, 'Alpha', 0.5);

    grid on;
    xlabel('频率 (GHz)', 'FontSize', 11);
    ylabel('\mu''_r', 'FontSize', 11);
    title('不同k值对应的相对磁导率实部', 'FontSize', 12);
    legend('Location', 'best', 'FontSize', 10);
    ylim([min(mu_r_real_all(:))-1, max(mu_r_real_all(:))+1]);

    %% 保存图片
    [~, filename, ~] = fileparts(s2p_file);
    fig_filename = sprintf('%s_diagnostic.png', filename);

    if options.verbose
        fprintf('\n保存图片: %s\n', fig_filename);
    end

    saveas(gcf, fig_filename);

    %% 保存Excel
    excel_filename = sprintf('%s_diagnostic.xlsx', filename);

    if options.verbose
        fprintf('保存Excel: %s\n', excel_filename);
    end

    % 准备数据表格
    % 列标题
    header = {'频率(GHz)', 'S11_幅度(dB)', 'S11_相位(deg)', ...
              'S21_幅度(dB)', 'S21_相位(deg)', ...
              '|Gamma|', 'Gamma相位(deg)', '|T|', 'T相位(deg)'};

    % 为每个k值添加εr和μr列
    for k_idx = 1:n_k
        k = k_values(k_idx);
        header{end+1} = sprintf('εr''(k=%+d)', k);
    end
    for k_idx = 1:n_k
        k = k_values(k_idx);
        header{end+1} = sprintf('μr''(k=%+d)', k);
    end

    % 数据
    data = [freq/1e9, ...
            20*log10(abs(S11)), angle(S11)*180/pi, ...
            20*log10(abs(S21)), angle(S21)*180/pi, ...
            abs(Gamma), angle(Gamma)*180/pi, ...
            abs(T), angle(T)*180/pi, ...
            epsilon_r_real_all, ...
            mu_r_real_all];

    % 写入Excel
    % 第一行：标题
    writecell(header, excel_filename, 'Sheet', 1, 'Range', 'A1');
    % 数据从第二行开始
    writematrix(data, excel_filename, 'Sheet', 1, 'Range', 'A2');

    if options.verbose
        fprintf('\n文件已保存:\n');
        fprintf('  图片: %s\n', fig_filename);
        fprintf('  数据: %s\n\n', excel_filename);
    end

    % 在命令行输出提示
    fprintf('请查看图中εr''和μr''的对比：\n');
    fprintf('  - 选择物理合理的曲线（平滑、符合材料特性）\n');
    fprintf('  - 空气：应选择εr''≈1且μr''≈1的k值\n');
    fprintf('  - 非磁性介电材料：应选择εr''>1且μr''≈1的k值\n');
    fprintf('  - 同时观察两个参数可以更准确地判断正确的k值\n\n');

end
