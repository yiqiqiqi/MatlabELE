function [S11, S21] = read_excel_s_params(excel_file, options)
% 从Excel读取S11和S21复数数据（不包含频率）
%
% 输入：
%   excel_file - Excel文件路径
%   options    - 可选参数
%       .sheet       - 工作表名称或编号，默认1
%       .start_row   - 数据起始行，默认1
%
% 输出：
%   S11  - S11复数向量
%   S21  - S21复数向量
%
% Excel格式（6列，从您的截图）：
%   列A: S11实部
%   列B: S11虚部（可能带'i'文本）
%   列C: （忽略）
%   列D: （忽略）
%   列E: S21实部
%   列F: S21虚部（可能带'i'文本）

    if nargin < 2
        options = struct();
    end

    if ~isfield(options, 'sheet')
        options.sheet = 1;
    end

    if ~isfield(options, 'start_row')
        options.start_row = 1;  % 默认从第1行开始（无标题）
    end

    fprintf('读取Excel文件: %s\n', excel_file);

    % 读取数值数据（自动跳过文本）
    try
        data = readmatrix(excel_file, 'Sheet', options.sheet);
    catch
        % 如果readmatrix失败，尝试xlsread
        try
            data = xlsread(excel_file, options.sheet);
        catch ME
            error('无法读取Excel文件: %s\n错误信息: %s', excel_file, ME.message);
        end
    end

    fprintf('  数据尺寸: %d行 × %d列\n', size(data, 1), size(data, 2));

    if size(data, 2) < 6
        error('Excel文件至少需要6列数据（S11实部、S11虚部、空、空、S21实部、S21虚部）');
    end

    % 从start_row开始读取
    data_rows = options.start_row:size(data, 1);
    n_points = length(data_rows);

    % 提取S11和S21
    % 列A: S11实部, 列B: S11虚部
    S11_real = data(data_rows, 1);
    S11_imag = data(data_rows, 2);

    % 列E: S21实部, 列F: S21虚部
    S21_real = data(data_rows, 5);
    S21_imag = data(data_rows, 6);

    % 组合成复数
    S11 = S11_real + 1j * S11_imag;
    S21 = S21_real + 1j * S21_imag;

    fprintf('  成功读取 %d 个数据点\n', n_points);
    fprintf('  S11范围: |S11| = %.4f ~ %.4f\n', min(abs(S11)), max(abs(S11)));
    fprintf('  S21范围: |S21| = %.4f ~ %.4f\n', min(abs(S21)), max(abs(S21)));

end
