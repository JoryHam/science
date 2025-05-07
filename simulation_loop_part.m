function [part_x, part_v, np, phi, w, global_x, history, history_count] = simulation_loop_part(...
    part_x, part_v, np, phi, w, w_s, ...
    nz, nr, ts, dh, Lz, Lr, spwt, q, v_drift, max_part, np_insert, dt, v_th, ...
    dist_type, beam_radius, beam_sigma_z, E_beam, B_type, start_it, ...
    history, history_count, snapshot_interval, total_steps)

    % 使用全局变量，避免重复定义
    global EPS0 A n0 phi0 Te M M_rel cell_volumes R_matrix efz efr B0 omega_c debug_mode 
    global Bz_grid Br_grid part_theta residual_shift use_single_precision c
    
    % 确保单精度一致性
    if use_single_precision
        if ~isa(part_x, 'single')
            part_x = single(part_x);
        end
        if ~isa(part_v, 'single')
            part_v = single(part_v);
        end
        if ~isa(w, 'single')
            w = single(w);
        end
        if ~isa(dh, 'single')
            dh = single(dh);
        end
        if ~isa(dt, 'single')
            dt = single(dt);
        end
        if ~isa(spwt, 'single')
            spwt = single(spwt);
        end
    end
    
    % 确保part_theta已初始化
    if ~exist('part_theta', 'var') || isempty(part_theta) || length(part_theta) < max_part
        fprintf('初始化粒子角度数组 part_theta (%d)\n', max_part);
        part_theta = zeros(max_part, 1);
        
        % 为现有粒子分配随机角度 (0到2π)
        if np > 0
            part_theta(1:np) = 2*pi*rand(np, 1);
        end
    end
    
    % 修改参数验证部分
    if ~exist('start_it', 'var')
        start_it = 1;  % 如果未提供，默认为1
    end
    
    % 初始化global_x数组
    global_x = zeros(max_part,1);
    if np > 0
        global_x(1:np) = part_x(1:np,1) + w;
    end
    
    % 确保时间步长合适
    dt_cycl = 2*pi/omega_c;
    dt = min(dt, dt_cycl/20);
    
    % 添加计时变量
    time_density = 0;
    time_field = 0;
    time_injection = 0;
    time_particle = 0;
    time_window = 0;
    time_boundary = 0;
    
    % 优化: 预计算常用值 - 只计算一次
    inv_dh = 1/dh;  % 避免循环中重复计算
    qm = q/M_rel;  % 使用相对论质量
    qm_dt_half = 0.5 * qm * dt;
    qm_dt = qm * dt;
    
    % 预分配数组 - 避免重复分配全局已存在的数组
    if use_single_precision
        E_z = zeros(max_part, 1, 'single');
        E_r = zeros(max_part, 1, 'single');
        Bz = zeros(max_part, 1, 'single');
        Br = zeros(max_part, 1, 'single');
        chg = zeros(nz, nr, 'single');
    else
        E_z = zeros(max_part, 1);
        E_r = zeros(max_part, 1);
        Bz = zeros(max_part, 1);
        Br = zeros(max_part, 1);
        chg = zeros(nz, nr);
    end
    
    % 在初始化部分添加或修改体积元素计算
    % 预计算体积元素
    if use_single_precision
        cell_volumes = zeros(nz, nr, 'single');
    else
        cell_volumes = zeros(nz, nr);
    end

    % 更精确的体积元素计算 - 使用环形体积公式
    for j = 1:nr
        % 对于r=0的特殊处理
        if j == 1
            r_outer = 0.5*dh;
            cell_volumes(:,j) = pi*r_outer^2*dh;
        else
            r_inner = (j-1.5)*dh;
            r_outer = (j-0.5)*dh;
            if r_inner < 0
                r_inner = 0;
            end
            % 使用环形体积公式
            cell_volumes(:,j) = pi*(r_outer^2 - r_inner^2)*dh;
        end
    end
    
    % 修改初始粒子位置检查
    if start_it == 1  % 只在第一次运行时执行
        invalid_particles = (part_x(:,1) < 0 | part_x(:,2) < 0);
        invalid_count = sum(invalid_particles);
        
        if invalid_count > 0
            if debug_mode
                fprintf('检测到%d个粒子位置无效，正在重新初始化...\n', invalid_count);
            end
            
            % 先将无效粒子移到数组末尾
            valid_np = np - invalid_count;
            if valid_np > 0
                valid_idx = find(~invalid_particles(1:np));
                
                % 保存有效粒子
                valid_part_x = part_x(valid_idx,:);
                valid_part_v = part_v(valid_idx,:);
                valid_part_theta = part_theta(valid_idx);
                valid_global_x = global_x(valid_idx);
                
                % 在数组前端放置有效粒子
                part_x(1:valid_np,:) = valid_part_x;
                part_v(1:valid_np,:) = valid_part_v;
                part_theta(1:valid_np) = valid_part_theta;
                global_x(1:valid_np) = valid_global_x;
                
                % 更新粒子计数
                np = valid_np;
            else
                % 如果所有粒子都无效，重置计数
                np = 0;
            end
            
            if debug_mode
                fprintf('有效粒子数: %d\n', np);
            end
        end
    end
    
    % 检查是否提供了history参数
    if nargin < 25
        % 如果未提供history参数，创建默认结构
        snapshot_interval = 50;  % 默认值
        history.part_x = cell(1, ceil(ts/snapshot_interval));
        history.part_v = cell(1, ceil(ts/snapshot_interval));
        history.part_theta = cell(1, ceil(ts/snapshot_interval));
        history.global_x = cell(1, ceil(ts/snapshot_interval));
        history.time_steps = zeros(1, ceil(ts/snapshot_interval));
        history_count = 0;
    end
    
    % 在仿真循环初始化部分修改
    if ~isfield(history, 'num_particles')
        history.num_particles = zeros(ts, 1);  % 使用ts参数（总时间步数）
    end

    if ~isfield(history, 'invalid_counts')
        history.invalid_counts = struct('boundary_z', zeros(ts, 1), ...
                                       'boundary_r', zeros(ts, 1), ...
                                       'invalid_data', zeros(ts, 1));
    end
    
    % 窗口位置使用双精度以提高长距离模拟的准确性
    if use_single_precision
        w = double(w);  % 窗口位置保持双精度
        w_s = double(w_s);  % 窗口移动量保持双精度
    end
    
    % 初始化并行环境（如果尚未初始化）
    if ~exist('pool_initialized', 'var') || isempty(pool_initialized)
        % 检查是否已有并行池
        if isempty(gcp('nocreate'))
            % 创建并行池，使用可用的核心数
            parpool('local', feature('numcores'));
        end
        pool_initialized = true;
    end
    
    % 主循环
    for it = 1:ts
        % 当前实际迭代步数（用于显示）
        current_it = it + start_it - 1;
        
        % 使用清零而不是重新分配
        chg(:) = 0;
        
        % 1. 计算电荷密度
        tic_local = tic;
        if np > 0
            % 向量化网格索引计算
            i = min(max(floor(part_x(1:np,1)*inv_dh + 0.5), 1), nz-1);
            j = min(max(floor(part_x(1:np,2)*inv_dh + 0.5), 1), nr-1);
            hx = part_x(1:np,1)*inv_dh - (i-0.5);
            hy = part_x(1:np,2)*inv_dh - (j-0.5);
            
            % 使用预先计算的权重
            for p = 1:np
                % 计算权重
                wp = (1-hx(p))*(1-hy(p));
                chg(i(p),j(p)) = chg(i(p),j(p)) + wp;
                
                wp = hx(p)*(1-hy(p));
                chg(i(p)+1,j(p)) = chg(i(p)+1,j(p)) + wp;
                
                wp = (1-hx(p))*hy(p);
                chg(i(p),j(p)+1) = chg(i(p),j(p)+1) + wp;
                
                wp = hx(p)*hy(p);
                chg(i(p)+1,j(p)+1) = chg(i(p)+1,j(p)+1) + wp;
            end
            
            % 转换为密度
            den = chg * spwt ./ cell_volumes;
        else
            den(:) = 0;
        end
        time_density = time_density + toc(tic_local);
        
        % 2. 求解电势
        tic_local = tic;
        A = setup_poisson_matrix(nz, nr, dh, R_matrix, w, Lz); % 动态更新矩阵
        phi = eval_2dpot_GS3(phi, w, Lz, dh);
        time_field = time_field + toc(tic_local);
        
        % 场插值开始 - 优化版本
        tic_interp = tic;
        
        % 1. 计算电场 - 使用向量化操作
        if ~exist('efz', 'var') || isempty(efz) || size(efz,1) ~= nz || size(efz,2) ~= nr
            % 只在需要时初始化电场数组
            if use_single_precision
                efz = zeros(nz, nr, 'single');
                efr = zeros(nz, nr, 'single');
            else
                efz = zeros(nz, nr);
                efr = zeros(nz, nr);
            end
        end
        
        % 向量化计算z方向电场 (内部点)
        efz(2:nz-1,:) = -(phi(3:nz,:) - phi(1:nz-2,:)) / (2*dh);
        
        % 边界处理
        efz(1,:) = -(phi(2,:) - phi(1,:)) / dh;
        efz(nz,:) = -(phi(nz,:) - phi(nz-1,:)) / dh;
        
        % 向量化计算r方向电场 (内部点)
        efr(:,2:nr-1) = -(phi(:,3:nr) - phi(:,1:nr-2)) / (2*dh);
        
        % 边界处理
        efr(:,1) = -(phi(:,2) - phi(:,1)) / dh;
        efr(:,nr) = -(phi(:,nr) - phi(:,nr-1)) / dh;
        
        % 2. 向量化场插值到粒子位置
        if np > 0
            % 计算网格索引
            iz = floor(part_x(1:np,1)/dh) + 1;
            ir = floor(part_x(1:np,2)/dh) + 1;
            
            % 边界检查
            iz = max(1, min(iz, nz-1));
            ir = max(1, min(ir, nr-1));
            
            % 计算插值权重
            wz = part_x(1:np,1)/dh - (iz-1);
            wr = part_x(1:np,2)/dh - (ir-1);
            
            % 预计算权重矩阵
            w00 = (1-wz).*(1-wr);
            w10 = wz.*(1-wr);
            w01 = (1-wz).*wr;
            w11 = wz.*wr;
            
            % 使用预计算的权重进行插值
            E_z(1:np) = zeros(np, 1);
            E_r(1:np) = zeros(np, 1);
            
            % 使用线性索引进行批量插值
            for i = 1:np
                E_z(i) = w00(i)*efz(iz(i),ir(i)) + w10(i)*efz(iz(i)+1,ir(i)) + ...
                         w01(i)*efz(iz(i),ir(i)+1) + w11(i)*efz(iz(i)+1,ir(i)+1);
                
                E_r(i) = w00(i)*efr(iz(i),ir(i)) + w10(i)*efr(iz(i)+1,ir(i)) + ...
                         w01(i)*efr(iz(i),ir(i)+1) + w11(i)*efr(iz(i)+1,ir(i)+1);
            end
            
            % 处理r=0轴上的特殊情况
            r_zero_mask = part_x(1:np,2) < dh/100;
            E_r(r_zero_mask) = 0;
            
            % 3. 如果有磁场，计算磁场
            if B0 > 0
                if strcmp(B_type, 'dipole')
                    % 计算偶极磁场
                    [Bz(1:np), Br(1:np)] = calculate_dipole_field(part_x(1:np,1), part_x(1:np,2), B0, Lz);
                else
                    % 均匀磁场
                    Bz(1:np) = B0;
                    Br(1:np) = 0;
                end
            end
        end
        
        interp_time = toc(tic_interp);
        if debug_mode
            fprintf('场插值时间: %.3f ms\n', interp_time*1000);
        end
        
        % 粒子推进开始
        tic_particle = tic;
        
        % 4. 向量化粒子推进
        if np > 0
            % 修改: 使用完整的Boris方法处理磁场
            if B0 > 0
                % 半步电场加速
                part_v(1:np,1) = part_v(1:np,1) + 0.5 * qm_dt * E_z(1:np);
                part_v(1:np,2) = part_v(1:np,2) + 0.5 * qm_dt * E_r(1:np);
                
                % 磁场旋转
                t_z = 0.5 * qm_dt * Bz(1:np);
                t_r = 0.5 * qm_dt * Br(1:np);
                t_theta = zeros(np, 1);  % 假设没有角向磁场
                
                % 计算t向量的模平方
                t_mag2 = t_z.^2 + t_r.^2;
                
                % 计算s向量
                s_z = 2 * t_z ./ (1 + t_mag2);
                s_r = 2 * t_r ./ (1 + t_mag2);
                s_theta = zeros(np, 1);  % 与t_theta对应
                
                % 保存原始速度
                v_minus_z = part_v(1:np,1);
                v_minus_r = part_v(1:np,2);
                v_minus_theta = part_v(1:np,3);
                
                % 计算v'
                v_prime_z = v_minus_z + v_minus_r .* t_theta - v_minus_theta .* t_r;
                v_prime_r = v_minus_r + v_minus_theta .* t_z - v_minus_z .* t_theta;
                v_prime_theta = v_minus_theta + v_minus_z .* t_r - v_minus_r .* t_z;
                
                % 计算v+
                part_v(1:np,1) = v_minus_z + v_prime_r .* s_theta - v_prime_theta .* s_r;
                part_v(1:np,2) = v_minus_r + v_prime_theta .* s_z - v_prime_z .* s_theta;
                part_v(1:np,3) = v_minus_theta + v_prime_z .* s_r - v_prime_r .* s_z;
                
                % 半步电场加速
                part_v(1:np,1) = part_v(1:np,1) + 0.5 * qm_dt * E_z(1:np);
                part_v(1:np,2) = part_v(1:np,2) + 0.5 * qm_dt * E_r(1:np);
            else
                % 无磁场情况下的简单推进
                part_v(1:np,1) = part_v(1:np,1) + qm_dt * E_z(1:np);
                part_v(1:np,2) = part_v(1:np,2) + qm_dt * E_r(1:np);
            end
            
            % 更新位置
            part_x(1:np,1) = part_x(1:np,1) + part_v(1:np,1) * dt;
            part_x(1:np,2) = part_x(1:np,2) + part_v(1:np,2) * dt;
            
            % 修改: 更新角向位置
            valid_r = part_x(1:np,2) > dh/100;  % 避免除以零
            if any(valid_r)
                part_theta(valid_r) = part_theta(valid_r) + part_v(valid_r,3) .* dt ./ part_x(valid_r,2);
                
                % 确保角度在0-2π范围内
                part_theta(1:np) = mod(part_theta(1:np), 2*pi);
            end
            
            % 修改: 处理粒子穿过r=0轴的情况
            r_negative = part_x(1:np,2) < 0;
            if any(r_negative)
                part_x(r_negative,2) = -part_x(r_negative,2);
                part_v(r_negative,2) = -part_v(r_negative,2);
                
                % 反射后角度增加π
                part_theta(r_negative) = mod(part_theta(r_negative) + pi, 2*pi);
            end
            
            % 更新全局位置
            global_x(1:np) = part_x(1:np,1) + w;
        end
        
        time_particle = time_particle + toc(tic_particle);
        if debug_mode
            fprintf('场插值时间: %.3f ms, 粒子推进时间: %.3f ms\n', interp_time*1000, time_particle*1000);
        end
        
        % 4.5 移动窗口 (在边界检查之前)
        tic_local = tic;
        if np > 0
            % 计算粒子在窗口中的平均位置
            mean_pos = mean(part_x(1:np,1));
            
            % 如果平均位置超过阈值，移动窗口
            if mean_pos > 0.7*Lz
                % 计算窗口移动量
                window_shift = 0.4*Lz;
                
                % 更新窗口位置
                w = w + window_shift;
                
                % 更新粒子局部坐标
                part_x(1:np,1) = part_x(1:np,1) - window_shift;
                
                % 更新全局坐标
                global_x(1:np) = part_x(1:np,1) + w;
                
                if debug_mode
                    fprintf('步数=%d: 移动窗口 %.2f m，新位置 %.2f m\n', current_it, window_shift, w);
                end

                % 移动窗口后检查z方向右边界
                tic_local = tic;
                for p = 1:np
                    % 只检查z方向右边界，左边界和r边界已在推进时检查
                    if part_x(p,1) >= Lz
                        % 标记粒子为无效
                        part_x(p,:) = -1;
                    end
                end
                time_boundary = time_boundary + toc(tic_local);
            end
        end
        time_window = toc(tic_local);
        
        % 5. 处理无效粒子
        tic_local = tic;
        if np > 0
            % 检查是否有NaN或无穷大
            invalid_data = isnan(part_x(1:np,1)) | isinf(part_x(1:np,1)) | ...
                          isnan(part_x(1:np,2)) | isinf(part_x(1:np,2)) | ...
                          isnan(part_v(1:np,1)) | isinf(part_v(1:np,1)) | ...
                          isnan(part_v(1:np,2)) | isinf(part_v(1:np,2)) | ...
                          isnan(part_v(1:np,3)) | isinf(part_v(1:np,3));
                          
            % 修改: 只检查z边界和r>Lr边界，不再检查r<0
            boundary_z = part_x(1:np,1) < 0 | part_x(1:np,1) > Lz;  % z边界
            boundary_r = part_x(1:np,2) > Lr;  % 只检查r>Lr边界
            
            % 合并所有无效标记
            invalid_idx = invalid_data | boundary_z | boundary_r;
            
            % 反转索引以找到有效粒子
            valid_idx = ~invalid_idx;
            
            % 统计无效粒子数
            n_invalid = sum(invalid_idx);
            
            % 输出信息
            if debug_mode && n_invalid > 0 && mod(it, 100) == 0
                fprintf('步数=%d: 移除了%d个无效粒子 (NaN=%d, Z边界=%d, R边界=%d)\n', ...
                       current_it, n_invalid, sum(invalid_data), ...
                       sum(boundary_z), sum(boundary_r));
            end
            
            % 更新粒子数组
            if n_invalid > 0
                if sum(valid_idx) > 0
                    part_x(1:sum(valid_idx),:) = part_x(valid_idx,:);
                    part_v(1:sum(valid_idx),:) = part_v(valid_idx,:);
                    part_theta(1:sum(valid_idx)) = part_theta(valid_idx);
                    global_x(1:sum(valid_idx)) = global_x(valid_idx);
                    np = sum(valid_idx);
                else
                    np = 0;
                end
            end
        end
        time_boundary = toc(tic_local);
        
        % 6. 注入新粒子
        tic_local = tic;
        if np_insert > 0 && np + np_insert <= max_part
            % 在z=0处注入粒子
            z_insert = 0.1 * dh * ones(np_insert, 1);
            
            % 根据分布类型选择不同的初始化方式
            switch dist_type
                case 'uniform'
                    % 真正的均匀分布（在XY平面上均匀）
                    r_insert = beam_radius * sqrt(rand(np_insert, 1));
                    % 计算角度
                    theta_insert = 2*pi*rand(np_insert, 1);
                    
                    % 为均匀分布设置速度
                    v_z = v_drift + v_th * randn(np_insert, 1);
                    v_r = v_th * randn(np_insert, 1);
                    v_theta = v_th * randn(np_insert, 1);
                    
                case 'gaussian'
                    % 空间分布保持，但使用反函数变换
                    r_insert = beam_radius * sqrt(-log(rand(np_insert, 1)));
                    theta_insert = 2*pi*rand(np_insert, 1);
                    
                    % 关键修改：不同方向使用不同的热速度
                    v_th_z = v_th;                       % 纵向热速度不变
                    v_th_r = v_th * 0.2;                 % 径向热速度降低为原来的20%
                    v_th_theta = v_th * 0.2;             % 角向热速度也降低
                    
                    % 纵向速度增加中心更高的梯度 (匹配文献中红色中心)
                    v_z = v_drift * (1 + 0.06*(1 - (r_insert/beam_radius).^2)) + v_th_z * randn(np_insert, 1);
                    
                    % 径向和角向速度使用缩小的热速度
                    v_r = v_th_r * randn(np_insert, 1);
                    v_theta = v_th_theta * randn(np_insert, 1);
                    
                otherwise
                    error('未知的分布类型: %s', dist_type);
            end
            
            % 相对论限制（对所有分布类型都应用）
            v_mag = sqrt(v_z.^2 + v_r.^2 + v_theta.^2);
            rel_limit = 0.999 * c;
            too_fast = v_mag > rel_limit;
            
            if any(too_fast)
                scale_factor = rel_limit ./ v_mag(too_fast);
                v_z(too_fast) = v_z(too_fast) .* scale_factor;
                v_r(too_fast) = v_r(too_fast) .* scale_factor;
                v_theta(too_fast) = v_theta(too_fast) .* scale_factor;
            end
            
            % 确保r>=0
            r_insert = max(r_insert, 0);
            
            % 添加新粒子到数组
            part_x(np+1:np+np_insert,:) = [z_insert, r_insert];
            part_v(np+1:np+np_insert,:) = [v_z, v_r, v_theta];
            part_theta(np+1:np+np_insert) = theta_insert;
            global_x(np+1:np+np_insert) = w + z_insert;
            
            % 更新粒子数量
            np = np + np_insert;
            
            if debug_mode && mod(it,100) == 0
                fprintf('步数=%d: 注入了%d个新粒子，总粒子数=%d\n', current_it, np_insert, np);
            end
        end
        time_injection = time_injection + toc(tic_local);
        
        % 7. 保存历史数据
        if mod(current_it, snapshot_interval) == 0 || it == ts
            % 记录当前窗口位置
            history_count = history_count + 1;
            history.z_positions(history_count) = w;
            
            % 计算束流RMS尺寸
            if np > 0
                % 计算束流中心
                r_mean = mean(part_x(1:np,2));
                % 计算RMS尺寸
                r_rms = sqrt(mean((part_x(1:np,2) - r_mean).^2));
                history.beam_rms_r(history_count) = r_rms;
                
                % 记录其他束流参数
                history.beam_mean_r(history_count) = r_mean;
                history.beam_max_r(history_count) = max(part_x(1:np,2));
                history.num_particles(history_count) = np;
                
                % 计算发散角(如果有足够的历史数据)
                if history_count > 1
                    dz = history.z_positions(history_count) - history.z_positions(history_count-1);
                    if dz > 0
                        dr = history.beam_rms_r(history_count) - history.beam_rms_r(history_count-1);
                        history.divergence_angle(history_count) = atan(dr/dz);
                    end
                end
            else
                % 如果没有粒子，记录NaN
                history.beam_rms_r(history_count) = NaN;
                history.beam_mean_r(history_count) = NaN;
                history.beam_max_r(history_count) = NaN;
                history.num_particles(history_count) = 0;
                history.divergence_angle(history_count) = NaN;
            end
            
            % 保存粒子数据
            history.part_x{history_count} = part_x(1:np,:);
            history.part_v{history_count} = part_v(1:np,:);
            history.part_theta{history_count} = part_theta(1:np);
            history.time_steps(history_count) = current_it;
            
            % 计算并保存束流RMS尺寸
            if np > 0
                beam_rms_r = std(part_x(1:np,2));  % 计算径向RMS尺寸
                history.beam_rms_r(history_count) = beam_rms_r;
                history.global_position(history_count) = w;  % 保存窗口位置
            else
                history.beam_rms_r(history_count) = 0;
                history.global_position(history_count) = w;
            end
            
            if debug_mode
                fprintf('保存历史数据点 #%d (步数=%d)，束流RMS尺寸=%.4e m，位置=%.2f m\n', ...
                        history_count, current_it, beam_rms_r, w);
            end
        end
        
        % 添加更详细的诊断
        if mod(it, 10) == 0 || np == 0
            % 计算径向统计量
            if np > 0
                r_mean = mean(part_x(1:np,2));
                r_max = max(part_x(1:np,2));
                r_std = std(part_x(1:np,2));
                fprintf('步数=%d: 粒子数=%d, 平均r=%.3e, 最大r=%.3e, 标准差=%.3e, 占边界比例=%.1f%%\n', ...
                       current_it, np, r_mean, r_max, r_std, 100*r_max/Lr);
            else
                fprintf('步数=%d: 警告 - 粒子数为零!\n', current_it);
            end
            
            % 如果粒子数为零，检查前一步的情况
            if np == 0 && exist('prev_part_x', 'var')
                r_prev = prev_part_x(:,2);
                fprintf('前一步统计: 粒子数=%d, 平均r=%.3e, 最大r=%.3e, r>Lr的比例=%.1f%%\n', ...
                       size(prev_part_x,1), mean(r_prev), max(r_prev), 100*sum(r_prev>Lr)/length(r_prev));
            end
            
            % 保存当前粒子位置用于下一步的诊断
            if np > 0
                prev_part_x = part_x(1:np,:);
            end
        end

        % 在每个时间步结束时添加
        history.num_particles(it) = np;

        % 在无效粒子处理代码中添加
        if exist('boundary_z', 'var') && (any(boundary_z) || any(boundary_r) || any(invalid_data))
            n_boundary_z = sum(boundary_z);
            n_boundary_r = sum(boundary_r);
            n_invalid_data = sum(invalid_data);
            
            % 记录不同类型的无效粒子数量
            history.invalid_counts.boundary_z(it) = n_boundary_z;
            history.invalid_counts.boundary_r(it) = n_boundary_r;
            history.invalid_counts.invalid_data(it) = n_invalid_data;
        end
    end
    
    % 输出性能统计
    if debug_mode
        total_time = time_density + time_field + time_injection + time_particle + time_window + time_boundary;
        fprintf('\n性能统计:\n');
        fprintf('密度计算: %.1f%% (%.3f s)\n', 100*time_density/total_time, time_density);
        fprintf('场计算: %.1f%% (%.3f s)\n', 100*time_field/total_time, time_field);
        fprintf('粒子注入: %.1f%% (%.3f s)\n', 100*time_injection/total_time, time_injection);
        fprintf('粒子推进: %.1f%% (%.3f s)\n', 100*time_particle/total_time, time_particle);
        fprintf('窗口移动: %.1f%% (%.3f s)\n', 100*time_window/total_time, time_window);
        fprintf('边界处理: %.1f%% (%.3f s)\n', 100*time_boundary/total_time, time_boundary);
        fprintf('平均每步耗时: %.3f ms\n', 1000*total_time/ts);
    end
    
    % 在仿真结束时保存最终结果
    if it == ts
        % 基本参数
        assignin('base', 'final_np', np);           % 最终粒子数
        assignin('base', 'final_w', w);             % 最终窗口位置
        
        % 物理量统计
        assignin('base', 'final_den_max', max(den(:)));    % 最终最大密度
        assignin('base', 'final_phi_max', max(phi(:)));    % 最终最大电势
        
        if np > 0
            v_mag = sqrt(sum(part_v(1:np,:).^2, 2));
            assignin('base', 'final_v_mean', mean(v_mag));  % 最终平均速度
            assignin('base', 'final_v_max', max(v_mag));    % 最终最大速度
        end
        
        % 在仿真结束时保存最终结果数据
        % 创建结构体存储最终状态
        final_state.den = den;
        final_state.phi = phi;
        final_state.E_mag = sqrt(efz.^2 + efr.^2);
        final_state.part_x = part_x(1:np,:);
        final_state.part_v = part_v(1:np,:);
        final_state.global_x = global_x(1:np);
        final_state.w = w;
        final_state.t = current_it*dt;
        
        % 保存到工作区
        assignin('base', 'final_state', final_state);
    end
end

% 只保留偶极磁场计算函数
function [Bz, Br] = calculate_dipole_field(z, r, B0, Lz)
    % 计算偶极磁场
    % z, r: 粒子位置
    % B0: 磁场强度参数
    % Lz: 计算域长度
    
    % 设置偶极子位置在计算域中心
    z0 = Lz/2;
    
    % 计算到偶极子的距离
    R = sqrt((z - z0).^2 + r.^2);
    % 防止除零
    R(R < 1e-10) = 1e-10;
    
    cos_theta = (z - z0) ./ R;  % cos(theta)
    sin_theta = r ./ R;         % sin(theta)
    
    % 偶极场分量 (B ~ 1/R^3)
    B_mag = B0 ./ (R.^3);
    Bz = B_mag .* (2 * cos_theta.^2 - sin_theta.^2);
    Br = B_mag .* (2 * cos_theta .* sin_theta);
    
    % 处理 R = 0 的情况，避免除以零
    Bz(R == 0) = B0;
    Br(R == 0) = 0;
end




