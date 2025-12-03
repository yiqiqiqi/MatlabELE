function [results] = calculate_from_s2p(s2p_file, d)
% 从S2P文件计算材料的电磁特性参数（X波段专用）
%
% 输入参数:
%   s2p_file: S2P文件路径
%   d: 样品厚度 (单位: 米)，例如 2e-3 表示 2mm
%
% 输出参数:
%   results: 结构体，包含所有计算结果
%       - frequency: 频率向量 (Hz)
%       - S11, S21: S参数（复数）
%       - Gamma: 反射系数
%       - T: 传输系数
%       - mu_r: 相对磁导率（复数）
%       - epsilon_r: 相对介电常数（复数）
%       - tan_delta: 损耗角正切
%
% 固定参数（X波段）:
%   - 截止波长 λc = 45.7 mm
%   - 波导宽边 a = 22.85 mm
%   - 适用频率范围: 8.2 - 12.4 GHz
%
% 使用示例:
%   results = calculate_from_s2p('material_test.s2p', 2e-3);
%
% 基于NRW方法 (方程 2-38 到 2-52)

    %% ========== 固定参数设置 ==========
    fprintf('========== X波段材料参数计算 ==========\n\n');

    % 物理常数
    c = 299792458;  % 光速 (m/s)
    mu_0 = 4*pi*1e-7;  % 真空磁导率 (H/m)
    epsilon_0 = 8.854187817e-12;  % 真空介电常数 (F/m)

    % X波段固定参数
    lambda_c = 45.7e-3;  % 截止波长 (m)
    waveguide_a = lambda_c / 2;  % 波导宽边 (m)

    % 假设入射介质为空气
    epsilon_i = 1;  % 相对介电常数
    mu_i = 1;       % 相对磁导率

    fprintf('固定参数:\n');
    fprintf('  截止波长 λc = %.2f mm\n', lambda_c * 1000);
    fprintf('  波导宽边 a = %.2f mm\n', waveguide_a * 1000);
    fprintf('  样品厚度 d = %.2f mm\n', d * 1000);
    fprintf('\n');

    %% ========== 读取S2P文件 ==========
    [freq, S11, S21, ~, ~] = read_s2p_file(s2p_file);

    n_points = length(freq);

    % 检查频率范围
    freq_ghz_min = min(freq) / 1e9;
    freq_ghz_max = max(freq) / 1e9;

    if freq_ghz_min < 8.0 || freq_ghz_max > 12.5
        warning('频率范围 (%.2f - %.2f GHz) 超出典型X波段范围 (8.2 - 12.4 GHz)', ...
                freq_ghz_min, freq_ghz_max);
    end

    %% ========== 计算电磁特性参数 ==========
    fprintf('开始计算材料参数...\n');

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

        % 如果两个都满足或都不满足，选择幅度较小的
        if (abs(Gamma_plus) <= 1) && (abs(Gamma_minus) <= 1)
            if abs(Gamma_minus) < abs(Gamma_plus)
                Gamma(i) = Gamma_minus;
            else
                Gamma(i) = Gamma_plus;
            end
        end

        % 方程 (2-42): 计算传输系数 T
        T(i) = (S11(i) + S21(i) - Gamma(i)) / (1 - (S11(i) + S21(i)) * Gamma(i));

        % 方程 (2-51): 计算 1/lambda^2
        % 需要处理对数的多值性
        log_T = log(1/T(i));

        % 如果虚部过大，可能需要调整（相位解缠）
        if abs(imag(log_T)) > pi
            % 简单的相位解缠
            log_T = log_T - 2*pi*1j*round(imag(log_T)/(2*pi));
        end

        inv_lambda_sq = -(1/(2*pi*d) * log_T)^2;
        lambda_g(i) = 1 / sqrt(inv_lambda_sq);

        % 方程 (2-48): 计算相对磁导率 mu_r
        term1 = sqrt(epsilon_i * mu_i / lambda_0(i)^2 - 1/lambda_c^2);

        % 检查是否在截止频率以下
        if imag(term1) ~= 0
            warning('频率点 %d (%.3f GHz) 可能接近或低于截止频率', i, freq(i)/1e9);
        end

        Lambda = d;  % 根据标准NRW方法

        numerator = (1 + Gamma(i)) * mu_i;
        denominator = Lambda * (1 - Gamma(i)) * term1;

        if abs(denominator) < 1e-10
            warning('分母接近零，频率点 %d', i);
            mu_r(i) = NaN + 1j*NaN;
        else
            mu_r(i) = numerator / denominator;
        end

        % 方程 (2-49): 计算相对介电常数 epsilon_r
        term2 = (inv_lambda_sq + 1/lambda_c^2) * lambda_0(i)^2;

        if abs(mu_r(i)) < 1e-10
            warning('μr接近零，频率点 %d', i);
            epsilon_r(i) = NaN + 1j*NaN;
        else
            epsilon_r(i) = term2 / mu_r(i);
        end

        % 方程 (2-50): 计算损耗角正切 tan(delta)
        if abs(real(epsilon_r(i))) > 1e-10
            tan_delta(i) = imag(epsilon_r(i)) / real(epsilon_r(i));
        else
            tan_delta(i) = NaN;
        end
    end

    fprintf('计算完成！\n\n');

    %% ========== 保存结果 ==========
    % 转换S参数为dB格式
    S11_dB = 20 * log10(abs(S11));
    S21_dB = 20 * log10(abs(S21));

    % 创建输出表格
    output_file = 'results_from_s2p.xlsx';

    output_table = table(freq/1e9, S11_dB, S21_dB, ...
                         abs(Gamma), angle(Gamma)*180/pi, ...
                         abs(T), angle(T)*180/pi, ...
                         real(mu_r), imag(mu_r), abs(mu_r), ...
                         real(epsilon_r), imag(epsilon_r), abs(epsilon_r), ...
                         tan_delta, ...
                         'VariableNames', {'Frequency_GHz', 'S11_dB', 'S21_dB', ...
                                          'Gamma_Magnitude', 'Gamma_Phase_deg', ...
                                          'T_Magnitude', 'T_Phase_deg', ...
                                          'mu_r_real', 'mu_r_imag', 'mu_r_magnitude', ...
                                          'epsilon_r_real', 'epsilon_r_imag', 'epsilon_r_magnitude', ...
                                          'tan_delta'});

    writetable(output_table, output_file);
    fprintf('结果已保存到: %s\n', output_file);

    %% ========== 显示统计信息 ==========
    fprintf('\n========== 计算结果统计 ==========\n\n');
    fprintf('频率范围: %.3f - %.3f GHz (%d 个点)\n', ...
            min(freq)/1e9, max(freq)/1e9, n_points);
    fprintf('\n相对磁导率 μr:\n');
    fprintf('  实部范围: %.4f ~ %.4f\n', min(real(mu_r)), max(real(mu_r)));
    fprintf('  虚部范围: %.4f ~ %.4f\n', min(imag(mu_r)), max(imag(mu_r)));
    fprintf('  平均值: %.4f + j%.4f\n', mean(real(mu_r)), mean(imag(mu_r)));

    fprintf('\n相对介电常数 εr:\n');
    fprintf('  实部范围: %.4f ~ %.4f\n', min(real(epsilon_r)), max(real(epsilon_r)));
    fprintf('  虚部范围: %.4f ~ %.4f\n', min(imag(epsilon_r)), max(imag(epsilon_r)));
    fprintf('  平均值: %.4f + j%.4f\n', mean(real(epsilon_r)), mean(imag(epsilon_r)));

    fprintf('\n损耗角正切 tanδ:\n');
    fprintf('  范围: %.6f ~ %.6f\n', min(tan_delta), max(tan_delta));
    fprintf('  平均值: %.6f\n', mean(tan_delta));

    %% ========== 绘制结果 ==========
    fprintf('\n正在绘制结果...\n');

    figure('Position', [100, 100, 1400, 900], 'Name', 'S2P材料参数分析结果');

    % 子图1: S参数
    subplot(3, 2, 1);
    plot(freq/1e9, S11_dB, 'b-', 'LineWidth', 1.5, 'DisplayName', 'S_{11}');
    hold on;
    plot(freq/1e9, S21_dB, 'r-', 'LineWidth', 1.5, 'DisplayName', 'S_{21}');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度 (dB)');
    legend('Location', 'best');
    title('S参数');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 子图2: 反射系数和传输系数
    subplot(3, 2, 2);
    plot(freq/1e9, abs(Gamma), 'b-', 'LineWidth', 1.5, 'DisplayName', '|\Gamma|');
    hold on;
    plot(freq/1e9, abs(T), 'r-', 'LineWidth', 1.5, 'DisplayName', '|T|');
    grid on;
    xlabel('频率 (GHz)');
    ylabel('幅度');
    legend('Location', 'best');
    title('反射系数和传输系数');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 子图3: 相对磁导率实部
    subplot(3, 2, 3);
    plot(freq/1e9, real(mu_r), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\mu_r''');
    title('相对磁导率（实部）');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 子图4: 相对磁导率虚部
    subplot(3, 2, 4);
    plot(freq/1e9, imag(mu_r), 'r-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\mu_r''''');
    title('相对磁导率（虚部）');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 子图5: 相对介电常数实部
    subplot(3, 2, 5);
    plot(freq/1e9, real(epsilon_r), 'b-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('\epsilon_r''');
    title('相对介电常数（实部）');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 子图6: 损耗角正切
    subplot(3, 2, 6);
    plot(freq/1e9, tan_delta, 'k-', 'LineWidth', 2);
    grid on;
    xlabel('频率 (GHz)');
    ylabel('tan\delta');
    title('损耗角正切');
    xlim([min(freq)/1e9, max(freq)/1e9]);

    % 保存图像
    saveas(gcf, 's2p_results_plot.png');
    fprintf('图像已保存到: s2p_results_plot.png\n');

    %% ========== 保存结果到结构体 ==========
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

    fprintf('\n========== 完成 ==========\n');

end
