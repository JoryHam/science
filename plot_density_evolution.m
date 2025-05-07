% 电子束密度演化可视化
% 绘制类似参考图像的8面板密度分布图

% 确保有历史数据
if ~exist('history', 'var') || ~exist('results', 'var')
    % 尝试加载最新的仿真结果
    files = dir('simulation_results_*.mat');
    if isempty(files)
        error('未找到仿真结果文件，请先运行仿真并保存结果');
    end
    
    % 按日期排序，找到最新的文件
    [~, idx] = sort([files.datenum], 'descend');
    latest_file = files(idx(1)).name;
    
    fprintf('加载最新的仿真结果：%s\n', latest_file);
    load(latest_file);
    
    if ~isfield(results, 'history')
        error('仿真结果中没有包含历史数据，请确保main.m保存了粒子历史');
    end
    
    history = results.history;
end

% 检查history结构
fprintf('检查history数据结构...\n');
if ~isfield(history, 'part_x')
    error('历史数据中缺少part_x字段');
end

% 检查part_x类型并显示信息
if iscell(history.part_x)
    fprintf('历史数据是元胞数组，包含%d个时间点\n', length(history.part_x));
else
    error('历史数据需要是元胞数组格式，每个元胞包含一个时间点的粒子位置');
end

% 从结果中提取物理参数，确保与main.m一致
if isfield(results, 'parameters')
    % 提取所有可用的物理参数
    params = results.parameters;
    
    % 显示主要物理参数
    fprintf('\n主要物理参数:\n');
    param_fields = fieldnames(params);
    for i = 1:length(param_fields)
        field = param_fields{i};
        if isnumeric(params.(field)) && isscalar(params.(field))
            fprintf('%s = %g\n', field, params.(field));
            % 将参数赋值到当前工作空间
            eval([field ' = params.' field ';']);
        end
    end
    
    % 特别处理一些关键参数
    if isfield(params, 'E_beam')
        E_beam = params.E_beam;
        fprintf('电子束能量: %.2f keV\n', E_beam/1e3);
    end
    
    if isfield(params, 'B0')
        B0 = params.B0;
        fprintf('磁场强度: %.3f T\n', B0);
    end
    
    % 计算相对论参数
    if isfield(params, 'E_beam')
        QE = 1.602e-19;
        M = 9.109e-31;
        c = 2.998e8;
        gamma = 1 + QE*E_beam/(M*c^2);
        v_drift = c*sqrt(1 - 1/gamma^2);
        fprintf('相对论因子 γ: %.3f\n', gamma);
        fprintf('漂移速度: %.2e m/s (%.2f%% c)\n', v_drift, 100*v_drift/c);
    end
else
    fprintf('警告: 结果中没有包含物理参数信息\n');
    % 使用默认值
    QE = 1.602e-19;
    M = 9.109e-31;
    c = 2.998e8;
    E_beam = 50e3;  % 默认50 keV
    gamma = 1 + QE*E_beam/(M*c^2);
    v_drift = c*sqrt(1 - 1/gamma^2);
end

% 设置固定的绘图范围（单位：米）
plot_range = 0.05;  % 50毫米 = 0.05米
plot_range_mm = plot_range * 1000;  % 转换为毫米，用于标尺显示

% 创建图像布局（2x4网格）
figure('Position', [100, 100, 1500, 700], 'Color', 'w');

% 选择8个时间点
ts = length(history.part_x);
fprintf('历史数据中共有%d个时间点\n', ts);

% 选择等间隔的时间点，并确保包括最后一个有效时间点
% 找到最后一个有效的时间点
last_valid = ts;
while last_valid > 0
    if ~isempty(history.part_x{last_valid}) && size(history.part_x{last_valid}, 2) >= 2
        break;
    end
    last_valid = last_valid - 1;
end

if last_valid < ts
    fprintf('注意：最后%d个时间点数据无效，使用时间点%d作为最后一个有效点\n', ts-last_valid, last_valid);
end

% 选择等间隔的时间点
time_indices = unique([round(linspace(1, last_valid-1, 7)) last_valid]);
if length(time_indices) < 8  % 如果有重复被移除
    time_indices = [time_indices ones(1, 8-length(time_indices))];
end
fprintf('选择的8个时间点索引：%s\n', mat2str(time_indices));

% 动态确定绘图范围（添加更多健壮性检查）
max_r = 0;
for i = 1:length(time_indices)
    t_idx = time_indices(i);
    if t_idx <= length(history.part_x)
        % 安全检查
        if ~isempty(history.part_x{t_idx}) && size(history.part_x{t_idx}, 2) >= 2
            cur_r = max(history.part_x{t_idx}(:,2));
            if ~isempty(cur_r) && ~isnan(cur_r)
                max_r = max(max_r, cur_r);
            end
        else
            fprintf('警告：时间点%d的数据为空或格式错误\n', t_idx);
        end
    else
        fprintf('警告：时间点索引%d超出范围\n', t_idx);
    end
end

% 绘制每个时间点的密度图
for i = 1:length(time_indices)
    t_idx = time_indices(i);
    subplot(2, 4, i);
    
    % 添加更多数据有效性检查
    if t_idx <= length(history.part_x) && ...
       ~isempty(history.part_x{t_idx}) && ...
       size(history.part_x{t_idx}, 2) >= 2 && ...
       t_idx <= length(history.part_v) && ...
       ~isempty(history.part_v{t_idx})
       
        % 提取数据
        part_x = history.part_x{t_idx};
        part_v = history.part_v{t_idx};
        np = size(part_x, 1);
        
        fprintf('时间点%d：处理%d个粒子\n', t_idx, np);
        
        if np > 0
            % 使用保存的粒子角度
            if isfield(history, 'part_theta') && ~isempty(history.part_theta{t_idx})
                theta = history.part_theta{t_idx};
            else
                % 如果没有保存的角度，使用之前的估计方法
                fprintf('未找到保存的角度数据，使用估计方法\n');
                
                % 使用速度方向估计
                vr = part_v(:,2);
                vtheta = part_v(:,3); 
                v_theta = atan2(vtheta, vr);
                
                % 添加时间相关的旋转和随机性
                t_idx_normalized = (t_idx-1) / (length(history.part_x)-1);
                if exist('omega_c', 'var') && ~isnan(omega_c) && omega_c ~= 0
                    omega = omega_c;
                else
                    QE = 1.602e-19;
                    M = 9.109e-31;
                    B0 = 0.5;
                    omega = QE*B0/M;
                end
                
                base_rotation = omega * t_idx_normalized * 20;
                theta = base_rotation + v_theta + 0.1*randn(size(v_theta));
            end
            
            % 转换为笛卡尔坐标
            r = part_x(:,2);
            x = r .* cos(theta);
            y = r .* sin(theta);
            
            % 修改：计算相对论动量比而不是速度比
            c = 2.998e8;  % 光速
            v_magnitudes = sqrt(sum(part_v.^2, 2));  % 速度大小
            
            % 计算相对论因子 gamma
            gamma = 1 ./ sqrt(1 - (v_magnitudes/c).^2);
            
            % 计算z方向动量 (p_z = gamma*m*v_z)
            pz = gamma .* part_v(:,1);
            
            % 使用从结果中提取的漂移速度计算参考动量
            gamma_0 = 1 / sqrt(1 - (v_drift/c)^2);
            p0 = gamma_0 * M * v_drift;
            fprintf('参考动量 p0: %.3e kg·m/s\n', p0);
            
            % 计算动量比
            ratio = pz * M / p0;  % 确保包含质量因子
            
            % 对大数据集进行采样
            max_points = 10000;
            if np > max_points
                idx = randperm(np, max_points);
                x = x(idx);
                y = y(idx);
                ratio = ratio(idx);
                fprintf('对时间点%d采样至%d个粒子\n', t_idx, max_points);
            end
            
            % 使用散点图带密度颜色
            scatter(x, y, 10, ratio, 'filled', 'MarkerFaceAlpha', 0.8);
            
            % 设置轴标签和网格
            xlabel('x (m)');
            ylabel('y (m)');
            grid on;
            set(gca, 'FontSize', 10);
            
            % 只在第一行第一列和第二行第一列子图显示y轴标签
            if ~(i == 1 || i == 5)
                set(gca, 'YTickLabel', []);
            end
            
            % 只在最后一行显示x轴标签
            if i <= 4
                set(gca, 'XTickLabel', []);
            end
            
            % 重要：确保所有子图使用完全相同的坐标轴范围
            % 移除 axis equal tight 命令，改为明确设置坐标轴范围
            axis equal;  % 保持纵横比例相等
            xlim([-plot_range, plot_range]);  % 强制X轴范围
            ylim([-plot_range, plot_range]);  % 强制Y轴范围
            
            % 可选：添加边框以便于区分子图
            box on;
        else
            text(0, 0, '无粒子数据', 'HorizontalAlignment', 'center');
        end
    else
        text(0, 0, '数据缺失或格式错误', 'HorizontalAlignment', 'center');
        fprintf('时间点%d的数据缺失或格式错误\n', t_idx);
    end
end

% 在循环外添加全局坐标轴范围检查
% 确保所有子图使用相同的坐标轴范围
all_axes = findall(gcf, 'type', 'axes');
for ax = 1:length(all_axes)
    if strcmp(get(all_axes(ax), 'Tag'), '')  % 跳过颜色条轴
        xlim(all_axes(ax), [-plot_range, plot_range]);
        ylim(all_axes(ax), [-plot_range, plot_range]);
    end
end

% 添加统一的色标
cb = colorbar('Position', [0.92, 0.1, 0.02, 0.8]);
ylabel(cb, 'p_z/p_0', 'FontSize', 12);

% 添加这行明确设置颜色轴范围
caxis([0.96, 1.04]);  % 与文献保持一致的物理意义范围

% 创建与文献匹配的自定义颜色映射
literature_colormap = [
    0.0  0.0  0.5;  % 0.96 - 深蓝
    0.0  0.0  0.8;  % 0.97 - 蓝色
    0.0  0.5  1.0;  % 0.98 - 浅蓝
    0.0  0.8  1.0;  % 0.99 - 淡蓝
    0.0  1.0  0.2;  % 1.00 - 绿色
    0.5  1.0  0.0;  % 1.01 - 淡黄绿
    1.0  1.0  0.0;  % 1.02 - 黄色
    1.0  0.5  0.0;  % 1.03 - 橙色
    1.0  0.0  0.0   % 1.04 - 红色
];

% 插值生成平滑的颜色映射
x = linspace(0, 1, size(literature_colormap, 1));
xq = linspace(0, 1, 256);
r = interp1(x, literature_colormap(:,1), xq, 'pchip');
g = interp1(x, literature_colormap(:,2), xq, 'pchip');
b = interp1(x, literature_colormap(:,3), xq, 'pchip');
smooth_colormap = [r' g' b'];

% 应用颜色映射
colormap(smooth_colormap);

% 添加总标题
sgtitle('电子束运动过程中横向尺寸和空间分布的演化情况', 'FontSize', 15);

% 添加标尺指示
annotation('textbox', [0.01, 0.01, 0.1, 0.05], 'String', ...
    sprintf('Scale: %d mm', plot_range*1000), 'EdgeColor', 'none');

% 添加物理参数标注到图像
param_text = sprintf('E_{beam}=%.0f keV, B=%.2f T', E_beam/1e3, B0);
annotation('textbox', [0.3, 0.01, 0.4, 0.03], 'String', param_text, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

% 保存图像选项
save_fig = input('是否保存密度演化图像? (y/n): ', 's');
if strcmpi(save_fig, 'y')
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('electron_beam_density_evolution_%s.png', timestamp);
    print(gcf, filename, '-dpng', '-r300');
    fprintf('密度演化图像已保存为: %s\n', filename);
end 