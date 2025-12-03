# 材料电磁特性参数计算工具

基于S参数（S11和S21）计算材料的相对磁导率、相对介电常数和损耗角正切。

## 功能特点

- ✅ **支持S2P文件**（Touchstone格式）- 从VNA直接导出 ⭐推荐
- ✅ 支持Excel文件读取S11和S21参数
- ✅ 计算反射系数Γ和传输系数T
- ✅ 计算相对磁导率μr（实部和虚部）
- ✅ 计算相对介电常数εr（实部和虚部）
- ✅ 计算损耗角正切tanδ
- ✅ 自动生成结果Excel文件和可视化图表

## 快速开始（S2P文件）⭐ 最简单！

如果你的数据是从**矢量网络分析仪（VNA）**直接导出的**.s2p文件**，这是最方便的方法：

```matlab
% 只需一行代码！
results = calculate_from_s2p('your_file.s2p', 2e-3);
%                              ↑ S2P文件      ↑ 样品厚度(米)
```

**固定参数（X波段）：**
- 截止波长 λc = 45.7 mm
- 波导宽边 a = 22.85 mm
- 适用频率：8.2 - 12.4 GHz

**就这么简单！** 程序会自动：
- 读取S2P文件中的频率和S参数
- 计算μr、εr、tanδ
- 生成Excel结果文件：`results_from_s2p.xlsx`
- 生成可视化图表：`s2p_results_plot.png`

详细示例请查看 `example_s2p_usage.m`

## 理论基础

本工具基于传输/反射法（NRW方法）的以下公式：

### 核心公式

**反射系数Γ** (方程 2-38 ~ 2-41):
```
X = (S₁₁² - S₂₁² + 1) / (2S₁₁)
Γ = X ± √(X² - 1)  (选择 |Γ| ≤ 1)
```

**传输系数T** (方程 2-42):
```
T = (S₁₁ + S₂₁ - Γ) / (1 - (S₁₁ + S₂₁)Γ)
```

**相对磁导率μr** (方程 2-48):
```
μr = (1 + Γ)μᵢ / [Λ(1 - Γ)√(εᵢμᵢ/λ₀² - 1/λc²)]
```

**相对介电常数εr** (方程 2-49):
```
εr = [(1/λ² + 1/λc²)λ₀²] / μr
其中: 1/λ² = -[1/(2πd) ln(1/T)]²
```

**损耗角正切** (方程 2-50):
```
tanδ = εr'' / εr'
```

## 文件说明

### S2P文件处理（⭐推荐用于VNA数据）

#### 1. `calculate_from_s2p.m` - **S2P文件专用计算函数**

**最方便的版本**，直接处理VNA导出的S2P文件。

**使用方法：**
```matlab
% 只需两个参数：文件路径 + 厚度
results = calculate_from_s2p('material.s2p', 2e-3);
```

**特点：**
- 自动读取频率信息（无需手动输入）
- 支持MA、DB、RI三种S2P格式
- 固定X波段参数（λc=45.7mm）
- 自动处理所有频率点

#### 2. `read_s2p_file.m` - **S2P文件读取函数**

底层文件读取函数，支持标准Touchstone格式。

#### 3. `example_s2p_usage.m` - **S2P使用示例**

完整的使用示例和说明文档。

---

### Excel文件处理（用于已处理的数据）

#### 4. `calculate_from_excel_simple.m`

**简单易用的Excel版本**，适合快速计算。

**使用方法：**
```matlab
% 在MATLAB命令窗口中运行：
calculate_from_excel_simple('your_data.xlsx')
```

**修改参数：**

打开文件，在"参数设置"部分修改以下参数：

```matlab
% 样品厚度 (单位: 米)
d = 2e-3;  % 例如: 2mm

% 波导宽边尺寸 (单位: 米)
waveguide_a = 22.86e-3;  % 例如: X波段 WR-90

% Excel列名
s11_column_name = 'S11_log';  % S11列名
s21_column_name = 'S21_log';  % S21列名
freq_column_name = 'SE_r';    % 频率列名
```

#### 5. `calculate_material_parameters.m`

**完整功能版本**，提供更多控制选项。

**使用方法：**
```matlab
% 设置参数
excel_file = 'your_data.xlsx';
d = 2e-3;  % 样品厚度 (m)
freq = 10e9;  % 频率 (Hz)
waveguide_params.a = 22.86e-3;  % 波导宽边 (m)
waveguide_params.mode = 'TE10';

% 调用函数
results = calculate_material_parameters(excel_file, d, freq, waveguide_params);

% 访问结果
disp(results.mu_r);       % 相对磁导率
disp(results.epsilon_r);  % 相对介电常数
disp(results.tan_delta);  % 损耗角正切
```

#### 6. `example_usage.m`

示例脚本，展示如何使用Excel完整功能版本。

## 快速开始

### 方法A：使用S2P文件（⭐推荐）

**步骤1：准备S2P文件**
- 从VNA（矢量网络分析仪）导出的标准.s2p文件
- 文件已包含频率信息，无需额外处理

**步骤2：运行计算**
```matlab
% 直接调用，只需设置样品厚度
results = calculate_from_s2p('your_file.s2p', 2e-3);
%                                                ↑
%                                           厚度(米): 2mm = 2e-3
```

**步骤3：查看结果**
- Excel文件：`results_from_s2p.xlsx`
- 图像文件：`s2p_results_plot.png`
- 命令窗口显示统计信息

**常用厚度参考：**
```matlab
0.5 mm  →  d = 0.5e-3
1.0 mm  →  d = 1.0e-3
2.0 mm  →  d = 2.0e-3
5.0 mm  →  d = 5.0e-3
10 mm   →  d = 10e-3
```

---

### 方法B：使用Excel文件

**步骤1：准备Excel文件**

确保Excel文件包含以下列：
- S11数据列（dB格式）
- S21数据列（dB格式）
- 频率列（可选，Hz）

示例Excel格式：
```
| Frequency  | S11_log | S21_log | ... |
|------------|---------|---------|-----|
| 8.0E+09    | -5.32   | -1.20   | ... |
| 8.0E+09    | -5.33   | -1.20   | ... |
| ...        | ...     | ...     | ... |
```

### 步骤2：修改参数

根据你的实验设置，修改以下参数：

1. **样品厚度 (d)**
   - 单位：米 (m)
   - 例如：2mm = 2e-3

2. **波导尺寸 (waveguide_a)**
   - 单位：米 (m)
   - 常用波导参数见下表

3. **Excel列名**
   - 根据实际Excel文件修改列名

### 步骤3：运行计算

```matlab
% 方法1：使用简化版本（推荐）
calculate_from_excel_simple('your_data.xlsx')

% 方法2：使用完整版本
results = calculate_material_parameters('your_data.xlsx', 2e-3, 10e9, struct('a', 22.86e-3));
```

### 步骤4：查看结果

- **Excel文件**: `calculated_results.xlsx` 或 `material_parameters_results.xlsx`
- **图像文件**: `results_plot.png` 或 `material_parameters_plot.png`
- **命令窗口**: 显示关键参数和计算结果

## 常用波导参数

| 波段 | 型号   | 宽边 a (mm) | 频率范围 (GHz) |
|------|--------|-------------|----------------|
| X    | WR-90  | 22.86       | 8.2 - 12.4     |
| Ku   | WR-62  | 15.80       | 12.4 - 18.0    |
| K    | WR-42  | 10.67       | 18.0 - 26.5    |
| Ka   | WR-28  | 7.11        | 26.5 - 40.0    |

## 输出结果

### Excel表格包含：
- 频率 (Frequency_Hz)
- S11和S21 (dB)
- 反射系数Γ的幅度和相位
- 传输系数T的幅度和相位
- 相对磁导率μr (实部、虚部、幅度)
- 相对介电常数εr (实部、虚部、幅度)
- 损耗角正切 tanδ

### 图表包含：
1. 相对磁导率 vs 频率/数据点
2. 相对介电常数 vs 频率/数据点
3. 损耗角正切 vs 频率/数据点
4. 反射系数和传输系数 vs 频率/数据点

## 注意事项

1. **单位转换**
   - 所有长度单位必须使用米 (m)
   - 频率单位必须使用赫兹 (Hz)
   - S参数假设为dB格式

2. **符号选择**
   - 程序自动选择 |Γ| ≤ 1 的解
   - 如果结果异常，请检查输入数据

3. **物理约束**
   - 反射系数 |Γ| ≤ 1
   - Re(1/T) > 0 （方程 2-52）

4. **Excel文件格式**
   - 确保Excel文件可以被MATLAB读取
   - 列名必须与代码中设置的一致
   - 数据不能有空值或非数字值

## 常见问题

### Q1: 找不到Excel列名怎么办？
运行代码时，程序会显示所有列名。根据显示的列名修改代码中的设置。

### Q2: 计算结果不合理怎么办？
检查以下项：
- 样品厚度d是否正确
- 波导尺寸是否正确
- S参数是否为dB格式
- 频率范围是否在波导工作范围内

### Q3: 如何处理多个频率点？
如果Excel中有频率列，程序会自动读取。否则需要在代码中设置频率向量。

### Q4: 相位是否需要解缠？
根据实际情况，如果相位跳变较大，可能需要相位解缠处理。

## 参考文献

本代码基于传输/反射法（Nicolson-Ross-Weir方法），具体理论参见：

1. A.M. Nicolson and G.F. Ross, "Measurement of the intrinsic properties of materials by time-domain techniques," IEEE Trans. Instrum. Meas., 1970.

2. W.B. Weir, "Automatic measurement of complex dielectric constant and permeability at microwave frequencies," Proc. IEEE, 1974.

## 技术支持

如有问题或建议，请联系开发者。

## 版本信息

- **版本**: 1.0
- **日期**: 2025-12-02
- **开发者**: Claude
- **MATLAB版本要求**: R2016a或更高版本

## 许可证

本代码仅供学习和研究使用。
