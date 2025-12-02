function [results] = calculate_material_parameters(excel_file, d, freq, waveguide_params)
% 计算材料的电磁特性参数（相对磁导率、相对介电常数、损耗角正切）
%
% 输入参数:
%   excel_file: Excel文件路径，包含S11和S21数据
%   d: 样品厚度 (m)
%   freq: 频率 (Hz)
%   waveguide_params: 波导参数结构体
%       - a: 波导宽边尺寸 (m)
%       - b: 波导窄边尺寸 (m) [可选]
%       - mode: 工作模式 (默认 'TE10')
%
% 输出参数:
%   results: 结构体，包含计算结果
%       - Gamma: 反射系数
%       - T: 传输系数
%       - mu_r: 相对磁导率 (复数)
%       - epsilon_r: 相对介电常数 (复数)
%       - tan_delta: 损耗角正切
%       - lambda: 波长
%       - lambda_g: 波导波长
%
% 基于方程 (2-38) 到 (2-52)

    % 物理常数
    c = 299792458;  % 光速 (m/s)
    mu_0 = 4*pi*1e-7;  % 真空磁导率 (H/m)
    epsilon_0 = 8.854187817e-12;  % 真空介电常数 (F/m)

    % 计算波长和截止波长
    lambda_0 = c ./ freq;  % 自由空间波长

    % 计算截止波长 (对于矩形波导TE10模式)
    if ~isfield(waveguide_params, 'mode')
        waveguide_params.mode = 'TE10';
    end

    if strcmp(waveguide_params.mode, 'TE10')
        lambda_c = 2 * waveguide_params.a;  % TE10模式截止波长
    else
        error('目前只支持TE10模式');
    end

    % 假设入射介质和测量介质为空气
    epsilon_i = 1;  % 相对介电常数
    mu_i = 1;       % 相对磁导率

    % 读取Excel文件中的S参数数据
    data = readtable(excel_file);

    % 提取S11和S21数据 (假设列名为 S11_log, S21_log 或类似)
    % 需要根据实际Excel文件的列名调整
    col_names = data.Properties.VariableNames;

    % 尝试找到S11和S21列
    s11_col = find_column(col_names, {'S11', 'S11_log'});
    s21_col = find_column(col_names, {'S21', 'S21_log'});

    if isempty(s11_col) || isempty(s21_col)
        error('无法在Excel文件中找到S11和S21列');
    end

    S11_log = data{:, s11_col};
    S21_log = data{:, s21_col};

    % 将对数形式转换为线性形式 (假设是dB)
    S11 = 10.^(S11_log / 20);
    S21 = 10.^(S21_log / 20);

    % 如果频率是向量，需要确保与S参数长度一致
    if length(freq) == 1
        freq = repmat(freq, length(S11), 1);
        lambda_0 = repmat(lambda_0, length(S11), 1);
    end

    % 初始化结果数组
    n_points = length(S11);
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    mu_r = zeros(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);
    lambda_g = zeros(n_points, 1);

    % 对每个频率点进行计算
    for i = 1:n_points
        % 方程 (2-40): 计算变量 X
        X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));

        % 方程 (2-41): 计算反射系数 Gamma
        % 需要选择符号，|Gamma| <= 1
        Gamma_temp = [X + sqrt(X^2 - 1), X - sqrt(X^2 - 1)];

        % 选择 |Gamma| <= 1 的解
        if abs(Gamma_temp(1)) <= 1
            Gamma(i) = Gamma_temp(1);
        else
            Gamma(i) = Gamma_temp(2);
        end

        % 方程 (2-51) 和 (2-42): 计算传输系数 T
        % T = (S11 + S21 - Gamma) / (1 - (S11 + S21) * Gamma)
        T(i) = (S11(i) + S21(i) - Gamma(i)) / (1 - (S11(i) + S21(i)) * Gamma(i));

        % 方程 (2-51): 计算 1/lambda^2
        inv_lambda_sq = -(1/(2*pi*d) * log(1/T(i)))^2;

        % 方程 (2-52): 判断符号是否正确
        if real(1/T(i)) <= 0
            warning('Re(1/T) <= 0 at point %d, 可能需要调整符号选择', i);
        end

        % 计算波导波长
        lambda_g(i) = 1 / sqrt(inv_lambda_sq);

        % 方程 (2-48): 计算相对磁导率 mu_r
        term1 = sqrt(epsilon_i * mu_i / lambda_0(i)^2 - 1/lambda_c^2);
        Lambda = d;  % 这里假设 Lambda = d，可能需要根据实际情况调整

        numerator = (1 + Gamma(i)) * mu_i;
        denominator = Lambda * (1 - Gamma(i)) * term1;
        mu_r(i) = numerator / denominator;

        % 方程 (2-49): 计算相对介电常数 epsilon_r
        term2 = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;
        epsilon_r(i) = term2 / mu_r(i);

        % 方程 (2-50): 计算损耗角正切 tan(delta)
        tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
    end

    % 保存结果到结构体
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
    results.frequency = freq;

    % 输出结果到Excel文件
    output_file = 'material_parameters_results.xlsx';
    output_table = table(freq, S11_log, S21_log, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(mu_r), imag(mu_r), abs(mu_r), ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, ...
                         'VariableNames', {'Frequency_Hz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Mag', 'Gamma_Phase_deg', ...
                                          'T_Mag', 'T_Phase_deg', ...
                                          'mu_r_real', 'mu_r_imag', 'mu_r_abs', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_abs', ...
                                          'tan_delta'});

    writetable(output_table, output_file);
    fprintf('计算完成！结果已保存到 %s\n', output_file);

    % 绘制结果图
    plot_results(results);
end

function col_idx = find_column(col_names, search_names)
    % 在列名中查找匹配的列
    col_idx = [];
    for i = 1:length(search_names)
        for j = 1:length(col_names)
            if contains(col_names{j}, search_names{i}, 'IgnoreCase', true)
                col_idx = j;
                return;
            end
        end
    end
end

function plot_results(results)
    % 绘制计算结果

    figure('Position', [100, 100, 1200, 800]);

    % 绘制相对磁导率
    subplot(2, 2, 1);
    plot(results.frequency/1e9, real(results.mu_r), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(results.frequency/1e9, imag(results.mu_r), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相对磁导率 \mu_r');
    legend('\mu_r^{\prime} (实部)', '\mu_r^{\prime\prime} (虚部)');
    title('相对磁导率 vs 频率');

    % 绘制相对介电常数
    subplot(2, 2, 2);
    plot(results.frequency/1e9, real(results.epsilon_r), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(results.frequency/1e9, imag(results.epsilon_r), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('相对介电常数 \epsilon_r');
    legend('\epsilon_r^{\prime} (实部)', '\epsilon_r^{\prime\prime} (虚部)');
    title('相对介电常数 vs 频率');

    % 绘制损耗角正切
    subplot(2, 2, 3);
    plot(results.frequency/1e9, results.tan_delta, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('损耗角正切 tan\delta');
    title('损耗角正切 vs 频率');

    % 绘制反射系数和传输系数
    subplot(2, 2, 4);
    plot(results.frequency/1e9, abs(results.Gamma), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(results.frequency/1e9, abs(results.T), 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('|\Gamma|', '|T|');
    title('反射系数和传输系数 vs 频率');

    % 保存图像
    saveas(gcf, 'material_parameters_plot.png');
    fprintf('图像已保存到 material_parameters_plot.png\n');
end
