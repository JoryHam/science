% 设置MATLAB性能优化选项
feature('JIT', 1);           % 启用JIT加速
maxNumCompThreads('auto');   % 自动设置线程数

% 首先清除所有变量和图形
clc;
clear all;

% 设置全局变量 - 整合所有需要在函数间共享的变量
global EPS0 QE den A n0 phi0 Te M M_rel cell_volumes R_matrix dh efz efr B0 omega_c debug_mode 
global Bz_grid Br_grid residual_shift part_theta use_single_precision c Z_mesh

% 添加单精度标志
use_single_precision = true;  % 设置为true启用单精度

% 设置调试模式
debug_mode = true;  % 设置调试模式

if debug_mode
    profile on;
    fprintf('初始化参数...\n');
end

% 物理常数 - 集中定义所有物理常数
EPS0 = 8.854e-12;        % 真空介电常数 (F/m)
QE = 1.602e-19;          % 基本电荷 (C)
M = 9.109e-31;           % 电子质量 (kg)
c = 2.998e8;             % 光速 (m/s)
k_B = 1.38e-23;          % 玻尔兹曼常数 (J/K)

% 输入设置
n0 = 1e16;               % 电子体密度 (m^-3)
Te = 1.0;                % 电子温度 (eV)
E_beam = 50e3;           % 电子束能量 (eV) = 50 keV
energy_spread = 0.01;    % 1% 能量分散

% 磁场设置
B_type = 'uniform';      % 'uniform' 或 'dipole'
B0 = 0;                  % 磁场强度 (T) = 0 Tesla

% 注入分布设置
dist_type = 'gaussian';  
beam_radius = 5e-3;      % 横向均方根尺寸 = 5 mm
beam_sigma_z = 5e-3;     % 纵向均方根尺寸 = 5 mm 
angle_rms = 0.001e-3;    % 发散角 = 0.001 mrad

% 在计算漂移速度时添加相对论修正
gamma = 1 + QE*E_beam/(M*c^2);  % 相对论因子 ≈ 1.098
v_drift = c*sqrt(1 - 1/gamma^2);  % 相对论修正后的速度 ≈ 1.23×10^8 m/s (41% c)
v_th = sqrt(2*QE*Te/M);          % 热速度 ≈ 5.93×10^5 m/s (0.48% 漂移速度)

% 检查初始速度是否合理
if v_drift >= c
    error('错误：初始漂移速度超过光速！');
end

if debug_mode
    fprintf('初始漂移速度: %.2e m/s (%.2f%% c)\n', v_drift, 100*v_drift/c);
    fprintf('相对论因子 γ: %.3f\n', gamma);
end

% 修正粒子质量
M_rel = M * gamma;  % 相对论质量 ≈ 1.00×10^-30 kg

% 在计算回旋频率时使用相对论质量
omega_c = QE*B0/M_rel;  % 修正后的回旋频率 = 0 (无磁场)

if debug_mode
    fprintf('电子温度: %.2f eV\n', Te);
    fprintf('计算得到的热速度: %.4e m/s (漂移速度的%.4f%%)\n', v_th, 100*v_th/v_drift);
end

% 设置模拟域 - 只定义一次网格参数
nz = 129;                % Z方向网格节点数
nr = 257;                % R方向网格节点数

% 计算德拜长度
Te_J = Te * QE;          % 转换为焦耳 ≈ 1.602×10^-19 J
lD = sqrt(EPS0*Te_J/(n0*QE^2));  % 德拜长度 ≈ 7.43×10^-5 m = 74.3 μm
dh = lD;                 % 网格间距 = 74.3 μm

% 计算域尺寸
Lz = (nz-1)*dh;          % Z方向长度 ≈ 9.51 mm
Lr = (nr-1)*dh;          % R方向长度 ≈ 19.1 mm

% 初始化网格 - 只初始化一次
[z_coords, r_coords] = deal(linspace(0,Lz,nz), linspace(0,Lr,nr));
[R_mesh, Z_mesh] = meshgrid(r_coords, z_coords);  % 注意顺序，使其符合圆柱坐标系约定

% 在调试信息中添加磁场类型
if debug_mode
    fprintf('磁场类型: %s\n', B_type);
    fprintf('磁场强度: %.2e T\n', B0);
end

% 计算各种物理时间尺度
plasma_freq = sqrt(n0*QE^2/(M*EPS0));  % 等离子体频率 ≈ 5.64×10^11 Hz

% 计算各种时间步长限制
dt_CFL = 0.5 * dh / v_drift;           % CFL条件 ≈ 3.02×10^-13 s
dt_plasma = 0.2 / plasma_freq;         % 等离子体振荡 ≈ 3.55×10^-13 s
dt_Debye = 0.2 * lD / v_th;            % 德拜时间 ≈ 2.50×10^-11 s

% 初始化时间步长数组
dt_array = [dt_CFL, dt_plasma, dt_Debye];

% 如果有磁场，添加回旋周期限制
if B0 > 0
    dt_cyclotron = 0.1 / omega_c;
    dt_array = [dt_array, dt_cyclotron];
end

% 选择最小的时间步长
dt = min(dt_array);  % 时间步长增大2倍

% 如果使用单精度，转换dt
if use_single_precision
    dt = single(dt);
end

% 第一阶段调用 - 使用固定窗口
target_distance = 20;  % 目标传播距离 = 20 米
total_steps = ceil(target_distance/(v_drift*dt));  % 总步数

% 估计运行时间
% 假设每步耗时3ms，总运行时间 ≈ 22.5小时

% 注入参数
np_insert = min((nr-1)*10, 1000);  % 限制每步注入的粒子数

% 计算通量和权重 - 修正为圆柱坐标系中的正确形式
% 在圆柱坐标系中，通量应该是通过圆面积的粒子数
flux = n0 * v_drift * pi * Lr^2;  % 入射粒子通量 (粒子/秒)
npt = flux * dt;                  % 每步实际粒子数
spwt = npt / np_insert;           % 宏粒子权重

% 计算预期总粒子数并设置最大粒子数
expected_particles = flux * (Lz/v_drift);
max_part = 50e2;  % 忽略预期计算，直接采用上限

if debug_mode
    fprintf('预期粒子数: %.2e\n', expected_particles);
    fprintf('最大粒子数: %.2e\n', max_part);
end

% 初始化粒子数组
if use_single_precision
    part_x = zeros(max_part, 2, 'single');  % [z,r]位置
    part_v = zeros(max_part, 3, 'single');  % [vz,vr,vθ]速度
    global_x = zeros(max_part, 1, 'single'); % 全局坐标
    part_theta = zeros(max_part, 1, 'single'); % 粒子角度
else
    part_x = zeros(max_part, 2);
    part_v = zeros(max_part, 3);
    global_x = zeros(max_part, 1);
    part_theta = zeros(max_part, 1);
end

% 添加这一行初始化np
np = 0;  % 初始粒子数为0

% 初始化场量 - 统一初始化，避免重复
if use_single_precision
    phi = zeros(nz, nr, 'single');
    den = zeros(nz, nr, 'single');
    efz = zeros(nz, nr, 'single');
    efr = zeros(nz, nr, 'single');
    Bz_grid = B0 * ones(nz, nr, 'single');  % 初始化为B0
    Br_grid = zeros(nz, nr, 'single');
    R_matrix = single(R_mesh);
else
    phi = zeros(nz, nr);
    den = zeros(nz, nr);
    efz = zeros(nz, nr);
    efr = zeros(nz, nr);
    Bz_grid = B0 * ones(nz, nr);  % 初始化为B0
    Br_grid = zeros(nz, nr);
    R_matrix = R_mesh;
end

% 计算单元体积 - 修正为圆柱坐标系中的正确形式
cell_volumes = calculate_cell_volumes_cylindrical(R_mesh, dh);

% 初始化系数矩阵A（用于求解泊松方程）
A = setup_poisson_matrix(nz, nr, dh, R_matrix, 0, Lz);

% 检查phi0是否被正确定义
phi0 = 0;                % 参考电势（标量）

if debug_mode
    fprintf('初始化完成。开始仿真...\n');
end

% 设置窗口参数
w = 0;                   % 初始窗口位置
w_s = v_drift*dt;        % 窗口速度等于粒子漂移速度

% 初始化残余移动量 - 用于精确跟踪窗口位置
residual_shift = 0;

% 在主循环前添加计时器
tic;

% 在主循环前添加性能设置
feature('JIT', 1);           % 确保JIT加速开启
feature('accel', 'on');      % 打开加速
set(0, 'DefaultFigureRenderer', 'opengl'); % 使用OpenGL加速图形

% ==========================================
% === 第一阶段：固定窗口仿真 ================
% ==========================================
ts_part1 = ceil(Lz / (v_drift*dt));  % 第一阶段步数：粒子穿过一个计算域的时间

% ==========================================
% === 第二阶段：移动窗口仿真 ================
% ==========================================
ts_part2 = total_steps - ts_part1;  % 剩余步数用于移动窗口阶段

% 在主循环开始前添加历史记录初始化 (确保在ts_part2定义之后)
snapshot_interval = max(floor(total_steps/50), 10);  % 确保不超过50个快照点
history.part_x = cell(1, ceil((ts_part1+ts_part2)/snapshot_interval));
history.part_v = cell(1, ceil((ts_part1+ts_part2)/snapshot_interval));
history.part_theta = cell(1, ceil((ts_part1+ts_part2)/snapshot_interval));  % 添加角度历史
history.time_steps = zeros(1, ceil((ts_part1+ts_part2)/snapshot_interval));
history.beam_rms_r = zeros(1, ceil((ts_part1+ts_part2)/snapshot_interval));  % 新增：束流RMS尺寸
history.global_position = zeros(1, ceil((ts_part1+ts_part2)/snapshot_interval));  % 新增：全局位置
history.z_positions = zeros(1, ceil(total_steps/snapshot_interval) + 1);
history.beam_mean_r = zeros(1, ceil(total_steps/snapshot_interval) + 1);
history.beam_max_r = zeros(1, ceil(total_steps/snapshot_interval) + 1);
history.num_particles = zeros(1, ceil(total_steps/snapshot_interval) + 1);
history.divergence_angle = zeros(1, ceil(total_steps/snapshot_interval) + 1);
history_count = 0;

if debug_mode
    fprintf('\n开始第一阶段：固定窗口仿真，共%d步\n', ts_part1);
end

% 在调试信息中添加分布类型
if debug_mode
    fprintf('束流分布类型: %s\n', dist_type);
    if strcmp(dist_type, 'gaussian')
        fprintf('高斯分布标准差: %.2e m\n', beam_sigma_z);
    else
        fprintf('均匀分布半径: %.2e m\n', beam_radius);
    end
end

% 在粒子注入前添加检查
if (np + np_insert) > max_part
    warning('粒子数达到上限 %d/%d', np, max_part);
    np_insert = max_part - np;
    if np_insert <= 0
        error('粒子存储已满，无法注入新粒子。请增加 max_part 或调整仿真参数。');
    end
    warning('注入数调整为 %d', np_insert);
end

% 在物理常数初始化后立即转换为单精度（移动到更早的位置）
if use_single_precision
    % 转换物理常数
    EPS0 = single(EPS0);
    QE = single(QE);
    M = single(M);
    c = single(c);
    
    % 转换模拟参数
    n0 = single(n0);
    Te = single(Te);
    E_beam = single(E_beam);
    energy_spread = single(energy_spread);
    B0 = single(B0);
    beam_radius = single(beam_radius);
    beam_sigma_z = single(beam_sigma_z);
    
    % 转换计算值
    v_drift = single(v_drift);
    gamma = single(gamma);
    omega_c = single(omega_c);
    v_th = single(v_th);
    lD = single(lD);
    dh = single(dh);
    dt = single(dt);
end

% 第一阶段调用 - 使用固定窗口
[part_x, part_v, np, phi, w, global_x, history, history_count] = simulation_loop_part(...
    part_x, part_v, np, phi, w, 0, ...  % 固定窗口，w_s=0
    nz, nr, ts_part1, dh, Lz, Lr, spwt, QE, ...
    v_drift, max_part, np_insert, dt, v_th, ...
    dist_type, beam_radius, beam_sigma_z, E_beam, B_type, 1, ...
    history, history_count, snapshot_interval);  % 添加历史参数

% 在第一阶段结束后显示信息
if debug_mode
    fprintf('第一阶段完成。开始移动窗口阶段...\n');
end

% 第二阶段调用 - 使用移动窗口 (w_s = v_drift*dt)
[part_x, part_v, np, phi, w, global_x, history, history_count] = simulation_loop_part(...
    part_x, part_v, np, phi, w, w_s, ...  % 移动窗口，w_s=v_drift*dt
    nz, nr, ts_part2, dh, Lz, Lr, spwt, QE, ...
    v_drift, max_part, np_insert, dt, v_th, ...
    dist_type, beam_radius, beam_sigma_z, E_beam, B_type, ts_part1+1, ...
    history, history_count, snapshot_interval);  % 添加历史参数

% 保存结果
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save_filename = sprintf('simulation_results_%s.mat', timestamp);

% 创建结果结构体
results = struct();
results.parameters = struct(...
    'n0', n0, ...
    'Te', Te, ...
    'E_beam', E_beam, ...
    'B0', B0, ...
    'dt', dt, ...
    'dh', dh, ...
    'grid_size', [nz, nr], ...
    'ts_part1', ts_part1, ...
    'ts_part2', ts_part2, ...
    'total_steps', total_steps);
results.final_state = struct(...
    'den', den, ...
    'phi', phi, ...
    'np', np);

% 如果配置为完整保存，则保存粒子数据
results.part_data = struct(...
    'part_x', part_x(1:np,:), ...
    'part_v', part_v(1:np,:), ...
    'global_x', global_x(1:np));

% 在保存结果部分添加历史数据 (results结构体创建后)
results.history = history;

save(save_filename, 'results', '-v7.3');  % 使用v7.3格式支持>2GB的变量

if debug_mode
    fprintf('结果已保存到: %s\n', save_filename);
end

% 显示仿真信息
if debug_mode
    profile viewer;
    sim_time = toc;
    fprintf('仿真完成。总耗时: %.2f秒\n', sim_time);
    fprintf('平均每步耗时: %.3f毫秒\n', sim_time*1000/total_steps);
    
    % 显示主要参数
    fprintf('\n主要参数:\n');
    fprintf('网格: %dx%d\n', nz, nr);
    fprintf('德拜长度: %.2e m\n', lD);
    fprintf('时间步长: %.2e s\n', dt);
    fprintf('磁场强度: %.2e T\n', B0);
    fprintf('粒子数: %d\n', np);
end

% 在初始化后添加物理参数检查
if debug_mode
    fprintf('\n物理参数检查：\n');
    fprintf('德拜长度: %.2e m\n', lD);
    fprintf('等离子体频率: %.2e Hz\n', plasma_freq);
    fprintf('回旋频率: %.2e Hz\n', omega_c);
    fprintf('回旋半径: %.2e m\n', v_drift/omega_c);
    fprintf('时间步长: %.2e s\n', dt);
end

% 在所有代码的最后添加清理代码
if debug_mode
    fprintf('\n仿真完全结束，清理资源...\n');
end

% 清理大数组以释放内存
clear part_x part_v global_x phi den efz efr;

% 添加结果分析和绘图
if isfield(history, 'beam_rms_r') && ~isempty(history.beam_rms_r)
    % 提取数据 - 确保只使用有效的历史记录
    z_positions = history.z_positions(1:history_count);
    beam_sizes = history.beam_rms_r(1:history_count);
    
    % 调试信息
    fprintf('总数据点: %d\n', length(z_positions));
    
    % 创建图像
    figure('Position', [100, 100, 600, 500]);
    
    % 绘制数据点
    plot(z_positions, beam_sizes, 'bo', 'MarkerSize', 6);
    hold on;
    
    % 线性拟合 - 使用中心化和缩放
    valid_idx = ~isnan(beam_sizes) & ~isnan(z_positions) & (z_positions >= 0);
    fprintf('有效数据点: %d\n', sum(valid_idx));
    
    if sum(valid_idx) > 2
        % 使用中心化和缩放
        z_valid = z_positions(valid_idx);
        b_valid = beam_sizes(valid_idx);
        
        % 计算均值和标准差
        z_mean = mean(z_valid);
        z_std = std(z_valid);
        if z_std == 0, z_std = 1; end
        
        fprintf('z均值: %.3e, z标准差: %.3e\n', z_mean, z_std);
        
        % 中心化和缩放
        z_scaled = (z_valid - z_mean) / z_std;
        
        % 使用缩放后的数据拟合
        [p, ~, mu] = polyfit(z_scaled, b_valid, 1);
        
        % 用于绘图的点
        z_fit = linspace(min(z_positions), max(z_positions), 100);
        z_fit_scaled = (z_fit - z_mean) / z_std;
        beam_fit = polyval(p, z_fit_scaled, [], mu);
        
        % 绘制拟合线
        plot(z_fit, beam_fit, 'b-', 'LineWidth', 1);
        
        % 计算束流发散角 (偏转角)
        if length(z_positions) > 1 && length(beam_sizes) > 1
            valid_idx = beam_sizes > 0;
            
            % 使用 MATLAB 内置的中心化和缩放选项
            [p, ~, mu] = polyfit(z_positions(valid_idx), beam_sizes(valid_idx), 1);
            
            % 计算斜率 (考虑缩放)
            slope = p(1);  % 已经是正确的斜率，因为 polyfit 内部处理了缩放
            beam_divergence = atan(slope);  % 弧度
            fprintf('束流发散角: %.4e rad (%.4f mrad)\n', beam_divergence, beam_divergence*1000);
            fprintf('z均值: %.3e, z标准差: %.3e\n', mu(1), mu(2));
            
            % 计算短距离发散角
            short_range_idx = z_positions <= 20 & beam_sizes > 0;
            if sum(short_range_idx) > 1
                [p_short, ~, mu_short] = polyfit(z_positions(short_range_idx), beam_sizes(short_range_idx), 1);
                beam_div_short = atan(p_short(1));  % 弧度
                fprintf('短距离束流发散角 (0-20m): %.4e rad (%.4f mrad)\n', beam_div_short, beam_div_short*1000);
            end
            
            % 计算长距离发散角 - 添加这部分代码
            long_range_idx = z_positions <= 100 & beam_sizes > 0;
            if sum(long_range_idx) > 1
                [p_long, ~, mu_long] = polyfit(z_positions(long_range_idx), beam_sizes(long_range_idx), 1);
                beam_div_long = atan(p_long(1));  % 弧度
                fprintf('长距离束流发散角 (0-100m): %.4e rad (%.4f mrad)\n', beam_div_long, beam_div_long*1000);
            end
            
            % 绘制拟合线
            z_fit = linspace(min(z_positions), max(z_positions), 100);
            z_fit_scaled = (z_fit - mu(1)) / mu(2);
            beam_fit = polyval(p, z_fit_scaled, [], mu);
            
            hold on;
            plot(z_fit, beam_fit, 'r-', 'LineWidth', 2);
            
            % 添加图例
            legend('模拟数据', sprintf('拟合线 (斜率 = %.4e)', slope), 'Location', 'NorthWest');
            
            % 在图上显示发散角
            text(0.1*max(z_positions), 0.8*max(beam_sizes), ...
                 sprintf('发散角 = %.4f mrad', beam_divergence*1000), ...
                 'FontSize', 12);
        end
    else
        fprintf('警告: 有效数据点不足，无法进行可靠的拟合\n');
    end
    
    % 设置图像属性
    grid on;
    xlabel('z/m', 'FontSize', 12);
    ylabel('\sigma/m', 'FontSize', 12);
    title('电子束尺寸随传播距离的变化', 'FontSize', 14);
    
    % 设置坐标轴范围
    if any(valid_idx)
        xlim([0, max(z_positions(valid_idx))]);
        ylim([0, max(1.2, max(beam_sizes(valid_idx))*1.2)]);
    end
    
    % 保存图像
    saveas(gcf, 'beam_size_vs_distance.png');
    saveas(gcf, 'beam_size_vs_distance.fig');
end

% 修改为圆柱坐标系中的正确体积计算函数
function cell_volumes = calculate_cell_volumes_cylindrical(R_mesh, dh)
    [nz, nr] = size(R_mesh);
    cell_volumes = zeros(nz, nr);
    
    % 在圆柱坐标系中，体积元素为 2*pi*r*dr*dz
    for i = 1:nz
        for j = 1:nr
            % 计算单元中心的r值
            r = (j-1)*dh;
            
            % 处理轴线上的特殊情况
            if j == 1
                % 轴线上的单元体积
                r_inner = 0;
                r_outer = 0.5*dh;
            else
                % 普通单元体积
                r_inner = r - 0.5*dh;
                r_outer = r + 0.5*dh;
                
                % 确保r_inner不为负
                r_inner = max(r_inner, 0);
            end
            
            % 计算z方向的因子
            z_factor = 1.0;
            if i == 1 || i == nz
                z_factor = 0.5;  % 边界上的单元只有一半体积
            end
            
            % 计算体积: 2*pi*r_avg*dr*dz
            % 使用环形体积公式: pi*(r_outer^2 - r_inner^2)*dz
            cell_volumes(i,j) = z_factor * pi * (r_outer^2 - r_inner^2) * dh;
        end
    end
end



