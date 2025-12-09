function [freq, S11, S21] = read_excel_s_params(excel_file, options)
% 从Excel读取S参数数据（复数格式）
%
% 输入：
%   excel_file - Excel文件路径
%   options    - 可选参数
%       .sheet       - 工作表名称或编号，默认1
%       .start_row   - 数据起始行，默认2（跳过标题）
%       .freq_unit   - 频率单位，'GHz'或'Hz'，默认'GHz'
%
% 输出：
%   freq - 频率向量（Hz）
%   S11  - S11复数向量
%   S21  - S21复数向量
%
% Excel格式假设（6列）：
%   列A: 频率
%   列B: S11实部
%   列C: S11虚部（可能带'i'或'j'文本）
%   列D: S21实部
%   列E: S21虚部（可能带'i'或'j'文本）
%   列F: （可能为空或其他数据）

    if nargin < 2
        options = struct();
    end

    if ~isfield(options, 'sheet')
        options.sheet = 1;
    end

    if ~isfield(options, 'start_row')
        options.start_row = 2;  % 默认第2行开始（跳过标题）
    end

    if ~isfield(options, 'freq_unit')
        options.freq_unit = 'GHz';
    end

    fprintf('读取Excel文件: %s\n', excel_file);

    % 读取所有数据（混合格式）
    try
        [~, ~, raw_data] = xlsread(excel_file, options.sheet);
    catch ME
        error('无法读取Excel文件: %s\n错误信息: %s', excel_file, ME.message);
    end

    % 确定数据范围
    n_rows = size(raw_data, 1);
    n_cols = size(raw_data, 2);

    fprintf('  数据尺寸: %d行 × %d列\n', n_rows, n_cols);

    if n_cols < 5
        error('Excel文件至少需要5列数据（频率、S11实部、S11虚部、S21实部、S21虚部）');
    end

    % 从start_row开始读取数据
    data_rows = options.start_row:n_rows;
    n_points = length(data_rows);

    freq = zeros(n_points, 1);
    S11 = zeros(n_points, 1);
    S21 = zeros(n_points, 1);

    fprintf('  开始解析数据...\n');

    for idx = 1:n_points
        row_idx = data_rows(idx);

        try
            % 读取频率（列A）
            freq_val = raw_data{row_idx, 1};
            if ischar(freq_val) || isstring(freq_val)
                freq(idx) = str2double(freq_val);
            else
                freq(idx) = freq_val;
            end

            % 读取S11（列B和C）
            S11_real = parse_value(raw_data{row_idx, 2});
            S11_imag = parse_value(raw_data{row_idx, 3});
            S11(idx) = S11_real + 1j * S11_imag;

            % 读取S21（列D和E，或E和F）
            % 先尝试列D和E
            if n_cols >= 5
                S21_real = parse_value(raw_data{row_idx, 4});
                S21_imag = parse_value(raw_data{row_idx, 5});
            else
                S21_real = 0;
                S21_imag = 0;
            end
            S21(idx) = S21_real + 1j * S21_imag;

        catch ME
            warning('第%d行数据解析失败: %s', row_idx, ME.message);
            % 使用前一个有效值或NaN
            if idx > 1
                freq(idx) = freq(idx-1);
                S11(idx) = S11(idx-1);
                S21(idx) = S21(idx-1);
            else
                freq(idx) = NaN;
                S11(idx) = NaN;
                S21(idx) = NaN;
            end
        end
    end

    % 频率单位转换
    if strcmpi(options.freq_unit, 'GHz')
        freq = freq * 1e9;  % GHz -> Hz
        fprintf('  频率范围: %.2f ~ %.2f GHz\n', min(freq)/1e9, max(freq)/1e9);
    else
        fprintf('  频率范围: %.2e ~ %.2e Hz\n', min(freq), max(freq));
    end

    fprintf('  成功读取 %d 个数据点\n', n_points);
    fprintf('  S11范围: |S11| = %.4f ~ %.4f\n', min(abs(S11)), max(abs(S11)));
    fprintf('  S21范围: |S21| = %.4f ~ %.4f\n', min(abs(S21)), max(abs(S21)));

end

function val = parse_value(cell_content)
    % 解析单元格内容，可能是数值、文本数值、或带'i'/'j'的复数文本

    if isnumeric(cell_content)
        val = cell_content;
        return;
    end

    if ischar(cell_content) || isstring(cell_content)
        str = char(cell_content);
        str = strtrim(str);  % 去除首尾空格

        % 移除末尾的'i'或'j'（如果存在）
        if endsWith(str, 'i') || endsWith(str, 'j')
            str = str(1:end-1);
            str = strtrim(str);
        end

        % 移除空格
        str = strrep(str, ' ', '');

        % 尝试转换为数值
        val = str2double(str);

        if isnan(val)
            warning('无法解析: "%s"，使用0', cell_content);
            val = 0;
        end
        return;
    end

    % 其他情况（例如空单元格）
    val = 0;
end
