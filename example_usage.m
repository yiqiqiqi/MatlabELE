%% 材料电磁特性参数计算 - 使用示例
% 本脚本演示如何使用 calculate_material_parameters 函数
% 基于S11和S21参数计算材料的相对磁导率、相对介电常数和损耗角正切

clear all;
close all;
clc;

%% 1. 设置输入参数

% Excel文件路径 (包含S11和S21数据)
% 注意：请根据实际文件路径修改
excel_file = 'your_data_file.xlsx';  % 修改为你的Excel文件路径

% 样品厚度 (单位: 米)
d = 2e-3;  % 例如: 2mm

% 测试频率 (单位: Hz)
% 可以是单个频率或频率向量
freq = 10e9;  % 例如: 10 GHz

% 如果Excel中有多个频率点，可以设置频率向量
% freq = linspace(8e9, 12e9, 100);  % 8-12 GHz，100个点

% 波导参数
waveguide_params.a = 22.86e-3;  % 波导宽边尺寸 (单位: 米)
                                % 例如: X波段 WR-90 波导 a = 22.86mm
waveguide_params.b = 10.16e-3;  % 波导窄边尺寸 (单位: 米) [可选]
waveguide_params.mode = 'TE10'; % 工作模式

%% 2. 调用计算函数
fprintf('开始计算材料电磁特性参数...\n');
fprintf('样品厚度: %.2f mm\n', d*1000);
fprintf('频率: %.2f GHz\n', freq/1e9);
fprintf('波导宽边: %.2f mm\n', waveguide_params.a*1000);
fprintf('截止波长: %.2f mm\n', 2*waveguide_params.a*1000);
fprintf('\n');

try
    results = calculate_material_parameters(excel_file, d, freq, waveguide_params);

    %% 3. 显示计算结果
    fprintf('\n========== 计算结果 ==========\n\n');

    % 显示第一个数据点的结果（如果有多个频率点）
    idx = 1;  % 可以修改为显示其他点

    fprintf('频率: %.3f GHz\n', results.frequency(idx)/1e9);
    fprintf('\n反射系数 Γ:\n');
    fprintf('  幅度: %.6f\n', abs(results.Gamma(idx)));
    fprintf('  相位: %.2f°\n', angle(results.Gamma(idx))*180/pi);

    fprintf('\n传输系数 T:\n');
    fprintf('  幅度: %.6f\n', abs(results.T(idx)));
    fprintf('  相位: %.2f°\n', angle(results.T(idx))*180/pi);

    fprintf('\n相对磁导率 μr:\n');
    fprintf('  实部 μr'': %.6f\n', real(results.mu_r(idx)));
    fprintf('  虚部 μr'''': %.6f\n', imag(results.mu_r(idx)));
    fprintf('  幅度 |μr|: %.6f\n', abs(results.mu_r(idx)));

    fprintf('\n相对介电常数 εr:\n');
    fprintf('  实部 εr'': %.6f\n', real(results.epsilon_r(idx)));
    fprintf('  虚部 εr'''': %.6f\n', imag(results.epsilon_r(idx)));
    fprintf('  幅度 |εr|: %.6f\n', abs(results.epsilon_r(idx)));

    fprintf('\n损耗角正切 tanδ: %.6f\n', results.tan_delta(idx));

    fprintf('\n波长信息:\n');
    fprintf('  自由空间波长 λ0: %.3f mm\n', results.lambda_0(idx)*1000);
    fprintf('  截止波长 λc: %.3f mm\n', results.lambda_c*1000);
    fprintf('  波导波长 λg: %.3f mm\n', abs(results.lambda_g(idx))*1000);

    fprintf('\n结果已保存到: material_parameters_results.xlsx\n');
    fprintf('图像已保存到: material_parameters_plot.png\n');

catch ME
    fprintf('错误: %s\n', ME.message);
    fprintf('请检查:\n');
    fprintf('  1. Excel文件路径是否正确\n');
    fprintf('  2. Excel文件中是否包含S11和S21列\n');
    fprintf('  3. 输入参数是否正确\n');
    rethrow(ME);
end

%% 4. 常用波导参数参考
fprintf('\n========== 常用波导参数参考 ==========\n');
fprintf('X波段 (WR-90):  a = 22.86 mm, 频率范围: 8.2-12.4 GHz\n');
fprintf('Ku波段 (WR-62): a = 15.80 mm, 频率范围: 12.4-18 GHz\n');
fprintf('K波段 (WR-42):  a = 10.67 mm, 频率范围: 18-26.5 GHz\n');
fprintf('Ka波段 (WR-28): a = 7.11 mm,  频率范围: 26.5-40 GHz\n');
