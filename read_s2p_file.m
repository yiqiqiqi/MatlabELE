function [freq, S11, S21, S12, S22] = read_s2p_file(filename)
% 读取S2P (Touchstone)格式文件
%
% 输入参数:
%   filename: S2P文件路径
%
% 输出参数:
%   freq: 频率向量 (Hz)
%   S11, S21, S12, S22: S参数（复数形式）
%
% 支持的格式:
%   - MA (幅度/角度)
%   - DB (dB/角度)
%   - RI (实部/虚部)

    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end

    % 初始化
    freq = [];
    S11 = [];
    S21 = [];
    S12 = [];
    S22 = [];

    % 默认参数
    freq_unit = 'GHz';  % 频率单位
    format_type = 'MA'; % 数据格式 (MA, DB, RI)

    fprintf('正在读取S2P文件: %s\n', filename);

    % 逐行读取
    line_num = 0;
    while ~feof(fid)
        line = fgetl(fid);
        line_num = line_num + 1;

        % 跳过空行
        if isempty(strtrim(line))
            continue;
        end

        % 处理选项行 (以 # 开头)
        if startsWith(strtrim(line), '#')
            tokens = strsplit(strtrim(line));

            % 解析频率单位
            for i = 1:length(tokens)
                token_upper = upper(tokens{i});
                if ismember(token_upper, {'HZ', 'KHZ', 'MHZ', 'GHZ'})
                    freq_unit = token_upper;
                end

                % 解析数据格式
                if ismember(token_upper, {'MA', 'DB', 'RI'})
                    format_type = token_upper;
                end
            end

            fprintf('  频率单位: %s\n', freq_unit);
            fprintf('  数据格式: %s\n', format_type);
            continue;
        end

        % 跳过注释行 (以 ! 开头)
        if startsWith(strtrim(line), '!')
            continue;
        end

        % 解析数据行
        try
            % 分割数据
            data = sscanf(line, '%f');

            % 检查数据点数量
            if length(data) < 9
                % 可能数据跨行，尝试读取下一行
                continue;
            end

            % 提取数据 (格式: freq S11_1 S11_2 S21_1 S21_2 S12_1 S12_2 S22_1 S22_2)
            f = data(1);
            s11_1 = data(2);
            s11_2 = data(3);
            s21_1 = data(4);
            s21_2 = data(5);
            s12_1 = data(6);
            s12_2 = data(7);
            s22_1 = data(8);
            s22_2 = data(9);

            % 转换频率单位到 Hz
            switch upper(freq_unit)
                case 'HZ'
                    f_hz = f;
                case 'KHZ'
                    f_hz = f * 1e3;
                case 'MHZ'
                    f_hz = f * 1e6;
                case 'GHZ'
                    f_hz = f * 1e9;
                otherwise
                    f_hz = f * 1e9; % 默认GHz
            end

            % 根据格式转换为复数
            switch upper(format_type)
                case 'MA'  % 幅度/角度（角度为度数）
                    s11_complex = s11_1 * exp(1j * s11_2 * pi/180);
                    s21_complex = s21_1 * exp(1j * s21_2 * pi/180);
                    s12_complex = s12_1 * exp(1j * s12_2 * pi/180);
                    s22_complex = s22_1 * exp(1j * s22_2 * pi/180);

                case 'DB'  % dB/角度
                    s11_complex = 10^(s11_1/20) * exp(1j * s11_2 * pi/180);
                    s21_complex = 10^(s21_1/20) * exp(1j * s21_2 * pi/180);
                    s12_complex = 10^(s12_1/20) * exp(1j * s12_2 * pi/180);
                    s22_complex = 10^(s22_1/20) * exp(1j * s22_2 * pi/180);

                case 'RI'  % 实部/虚部
                    s11_complex = s11_1 + 1j * s11_2;
                    s21_complex = s21_1 + 1j * s21_2;
                    s12_complex = s12_1 + 1j * s12_2;
                    s22_complex = s22_1 + 1j * s22_2;

                otherwise
                    error('未知的数据格式: %s', format_type);
            end

            % 添加到数组
            freq = [freq; f_hz];
            S11 = [S11; s11_complex];
            S21 = [S21; s21_complex];
            S12 = [S12; s12_complex];
            S22 = [S22; s22_complex];

        catch ME
            % 忽略无法解析的行
            continue;
        end
    end

    % 关闭文件
    fclose(fid);

    % 检查是否读取到数据
    if isempty(freq)
        error('未能从文件中读取到有效数据');
    end

    fprintf('  成功读取 %d 个频率点\n', length(freq));
    fprintf('  频率范围: %.3f - %.3f GHz\n', min(freq)/1e9, max(freq)/1e9);
    fprintf('\n');

end
