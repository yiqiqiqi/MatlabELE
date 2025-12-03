function calculate_from_excel_simple(excel_file)
% 简化版本：直接从Excel文件计算材料电磁特性参数
%
% 使用方法:
%   calculate_from_excel_simple('your_file.xlsx')
%
% 输入参数:
%   excel_file: 包含S11和S21数据的Excel文件路径
%
% 注意：请在代码中修改以下参数以匹配你的实验设置:
%   - d: 样品厚度
%   - waveguide_a: 波导宽边尺寸
%   - S11和S21的列名

    %% ========== 参数设置 ==========
    % 请根据实际情况修改以下参数

    % 样品厚度 (单位: 米)
    d = 2e-3;  % 默认: 2mm，请修改为实际厚度

    % 波导宽边尺寸 (单位: 米)
    waveguide_a = 22.86e-3;  % 默认: X波段 WR-90 (22.86mm)

    % S参数列名（请根据Excel文件中的实际列名修改）
    s11_column_name = 'S11_log';  % S11列名
    s21_column_name = 'S21_log';  % S21列名
    freq_column_name = 'SE_r';    % 频率列名（如果有）

    % 如果Excel中没有频率列，请手动设置频率
    manual_freq = 10e9;  % 单位: Hz，例如 10 GHz
    use_manual_freq = false;  % 如果Excel中有频率数据，设置为false

    %% ========== 物理常数 ==========
    c = 299792458;  % 光速 (m/s)
    mu_0 = 4*pi*1e-7;  % 真空磁导率 (H/m)
    epsilon_0 = 8.854187817e-12;  % 真空介电常数 (F/m)

    % 假设入射介质和测量介质为空气
    epsilon_i = 1;  % 相对介电常数
    mu_i = 1;       % 相对磁导率

    % 计算截止波长 (TE10模式)
    lambda_c = 2 * waveguide_a;

    %% ========== 读取Excel数据 ==========
    fprintf('正在读取Excel文件: %s\n', excel_file);
    data = readtable(excel_file);

    % 显示列名，帮助用户确认
    fprintf('\nExcel文件中的列名:\n');
    disp(data.Properties.VariableNames);

    % 提取S参数数据
    try
        % 尝试读取S11和S21 (假设是dB格式)
        if ismember(s11_column_name, data.Properties.VariableNames)
            S11_dB = data.(s11_column_name);
        else
            error('找不到列: %s', s11_column_name);
        end

        if ismember(s21_column_name, data.Properties.VariableNames)
            S21_dB = data.(s21_column_name);
        else
            error('找不到列: %s', s21_column_name);
        end

        % 转换为线性值
        S11 = 10.^(S11_dB / 20);
        S21 = 10.^(S21_dB / 20);

        % 读取频率数据
        if use_manual_freq
            freq = repmat(manual_freq, length(S11), 1);
        else
            if ismember(freq_column_name, data.Properties.VariableNames)
                freq = data.(freq_column_name);
            else
                warning('找不到频率列: %s，使用手动设置的频率', freq_column_name);
                freq = repmat(manual_freq, length(S11), 1);
            end
        end

    catch ME
        fprintf('\n错误: %s\n', ME.message);
        fprintf('\n请检查并修改代码中的列名设置:\n');
        fprintf('  s11_column_name = ''%s'';\n', s11_column_name);
        fprintf('  s21_column_name = ''%s'';\n', s21_column_name);
        fprintf('  freq_column_name = ''%s'';\n', freq_column_name);
        rethrow(ME);
    end

    n_points = length(S11);
    fprintf('\n读取到 %d 个数据点\n', n_points);

    %% ========== 计算电磁特性参数 ==========
    fprintf('\n开始计算...\n');

    % 计算自由空间波长
    lambda_0 = c ./ freq;

    % 初始化结果数组
    Gamma = zeros(n_points, 1);
    T = zeros(n_points, 1);
    mu_r = zeros(n_points, 1);
    epsilon_r = zeros(n_points, 1);
    tan_delta = zeros(n_points, 1);
    lambda_g = zeros(n_points, 1);

    % 逐点计算
    for i = 1:n_points
        % 方程 (2-40): 计算变量 X
        X = (S11(i)^2 - S21(i)^2 + 1) / (2 * S11(i));

        % 方程 (2-41): 计算反射系数 Gamma
        % 选择 |Gamma| <= 1 的解
        Gamma_plus = X + sqrt(X^2 - 1);
        Gamma_minus = X - sqrt(X^2 - 1);

        if abs(Gamma_plus) <= 1
            Gamma(i) = Gamma_plus;
        else
            Gamma(i) = Gamma_minus;
        end

        % 方程 (2-42): 计算传输系数 T
        T(i) = (S11(i) + S21(i) - Gamma(i)) / (1 - (S11(i) + S21(i)) * Gamma(i));

        % 方程 (2-51): 计算 1/lambda^2
        inv_lambda_sq = -(1/(2*pi*d) * log(1/T(i)))^2;
        lambda_g(i) = 1 / sqrt(inv_lambda_sq);

        % 方程 (2-48): 计算相对磁导率 mu_r
        term1 = sqrt(epsilon_i * mu_i / lambda_0(i)^2 - 1/lambda_c^2);
        Lambda = d;  % 根据实际情况，Lambda通常等于样品厚度d

        numerator = (1 + Gamma(i)) * mu_i;
        denominator = Lambda * (1 - Gamma(i)) * term1;
        mu_r(i) = numerator / denominator;

        % 方程 (2-49): 计算相对介电常数 epsilon_r
        term2 = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;
        epsilon_r(i) = term2 / mu_r(i);

        % 方程 (2-50): 计算损耗角正切 tan(delta)
        if real(epsilon_r(i)) ~= 0
            tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
        else
            tan_delta(i) = NaN;
        end
    end

    fprintf('计算完成！\n');

    %% ========== 保存结果 ==========
    output_file = 'calculated_results.xlsx';

    output_table = table(freq, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(mu_r), imag(mu_r), abs(mu_r), ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, ...
                         'VariableNames', {'Frequency_Hz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Magnitude', 'Gamma_Phase_deg', ...
                                          'T_Magnitude', 'T_Phase_deg', ...
                                          'mu_r_real', 'mu_r_imag', 'mu_r_magnitude', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_magnitude', ...
                                          'tan_delta'});

    writetable(output_table, output_file);
    fprintf('\n结果已保存到: %s\n', output_file);

    %% ========== 显示部分结果 ==========
    fprintf('\n========== 计算结果示例（前5个点） ==========\n');
    for i = 1:min(5, n_points)
        fprintf('\n--- 数据点 %d ---\n', i);
        if ~use_manual_freq
            fprintf('频率: %.3f GHz\n', freq(i)/1e9);
        end
        fprintf('S11: %.2f dB, S21: %.2f dB\n', S11_dB(i), S21_dB(i));
        fprintf('相对磁导率 μr = %.4f + j%.4f (|μr| = %.4f)\n', ...
                real(mu_r(i)), imag(mu_r(i)), abs(mu_r(i)));
        fprintf('相对介电常数 εr = %.4f + j%.4f (|εr| = %.4f)\n', ...
                real(epsilon_r(i)), imag(epsilon_r(i)), abs(epsilon_r(i)));
        fprintf('损耗角正切 tanδ = %.6f\n', tan_delta(i));
    end

    %% ========== 绘制结果 ==========
    fprintf('\n正在绘制结果...\n');

    figure('Position', [100, 100, 1200, 800], 'Name', '材料电磁特性参数');

    % 绘制相对磁导率
    subplot(2, 2, 1);
    if use_manual_freq
        x_data = 1:n_points;
        x_label = '数据点';
    else
        x_data = freq/1e9;
        x_label = '频率 (GHz)';
    end

    plot(x_data, real(mu_r), 'b-', 'LineWidth', 1.5, 'DisplayName', '\mu_r'' (实部)');
    hold on;
    plot(x_data, imag(mu_r), 'r--', 'LineWidth', 1.5, 'DisplayName', '\mu_r'''' (虚部)');
    grid on;
    xlabel(x_label);
    ylabel('相对磁导率 \mu_r');
    legend('Location', 'best');
    title('相对磁导率');

    % 绘制相对介电常数
    subplot(2, 2, 2);
    plot(x_data, real(epsilon_r), 'b-', 'LineWidth', 1.5, 'DisplayName', '\epsilon_r'' (实部)');
    hold on;
    plot(x_data, imag(epsilon_r), 'r--', 'LineWidth', 1.5, 'DisplayName', '\epsilon_r'''' (虚部)');
    grid on;
    xlabel(x_label);
    ylabel('相对介电常数 \epsilon_r');
    legend('Location', 'best');
    title('相对介电常数');

    % 绘制损耗角正切
    subplot(2, 2, 3);
    plot(x_data, tan_delta, 'k-', 'LineWidth', 1.5);
    grid on;
    xlabel(x_label);
    ylabel('损耗角正切 tan\delta');
    title('损耗角正切');

    % 绘制反射系数和传输系数
    subplot(2, 2, 4);
    plot(x_data, abs(Gamma), 'b-', 'LineWidth', 1.5, 'DisplayName', '|\Gamma|');
    hold on;
    plot(x_data, abs(T), 'r-', 'LineWidth', 1.5, 'DisplayName', '|T|');
    grid on;
    xlabel(x_label);
    ylabel('幅度');
    legend('Location', 'best');
    title('反射系数和传输系数');

    % 保存图像
    saveas(gcf, 'results_plot.png');
    fprintf('图像已保存到: results_plot.png\n');

    fprintf('\n========== 完成 ==========\n');
    fprintf('\n参数设置:\n');
    fprintf('  样品厚度 d = %.2f mm\n', d*1000);
    fprintf('  波导宽边 a = %.2f mm\n', waveguide_a*1000);
    fprintf('  截止波长 λc = %.2f mm\n', lambda_c*1000);
end
