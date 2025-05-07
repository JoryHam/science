function [x] = eval_2dpot_GS3(phi, w, Lz, dh)
    tic_local = tic;
    % 声明全局变量
    global EPS0 QE den A n0 phi0 Te M cell_volumes R_matrix debug_mode use_single_precision
    
    % 获取网格尺寸
    [nz, nr] = size(den);
    nn = nz * nr;
    
    % 使用上一步的解作为初始猜测
    x = phi;
    
    % 计算系数
    coef = -QE/EPS0;
    
    % 计算右端项 (使用预计算的体积元素)
    b = den(:) .* cell_volumes(:) * coef;
    
    % 应用边界条件到右端项 (向量化操作)
    % 只对角点应用边界条件
    b(1) = phi0;  % 左下角 (z=0, r=0)
    b(nz) = 0;    % 右下角 (z=Lz, r=0)
    
    % 外边界 (r=Lr)
    b((nr-1)*nz+1:end) = 0;
    
    % z边界
    for j = 1:nr
        b((j-1)*nz+1) = phi0;  % z=0
        b(j*nz) = 0;           % z=Lz
    end
    
    % 确保b是双精度，因为稀疏矩阵运算不支持单精度
    if use_single_precision
        b = double(b);
    end
    
    % 使用缓存的系数矩阵和LU分解
    persistent A_cached L U P Q last_size;
    
    % 只在第一次或网格尺寸变化时重建矩阵和分解
    if isempty(A_cached) || isempty(last_size) || last_size ~= nn
        if debug_mode
            fprintf('重建泊松矩阵和LU分解...\n');
        end
        
        % 复制系数矩阵
        A_cached = A;
        
        % 应用边界条件
        % z=0 边界
        for j = 1:nr
            u = (j-1)*nz + 1;
            A_cached(u,:) = 0;
            A_cached(u,u) = 1;
        end
        
        % z=Lz 边界
        for j = 1:nr
            u = (j-1)*nz + nz;
            A_cached(u,:) = 0;
            A_cached(u,u) = 1;
        end
        
        % 轴对称边界 r=0 (修正物理正确性)
        for i = 2:nz-1
            u = i;  % 轴上的点，j=1
            A_cached(u,:) = 0;
            % 在r=0处，使用L'Hôpital法则处理1/r项
            A_cached(u,u) = -4/dh^2 - 2/dh^2;  % 修正中心点系数
            A_cached(u,u+nz) = 4/dh^2;         % r+1点，系数为4/dh^2
            A_cached(u,u-1) = 1/dh^2;          % z-1
            A_cached(u,u+1) = 1/dh^2;          % z+1
        end
        
        % 处理轴上的角点
        u = 1;  % 左下角 (z=0, r=0)
        A_cached(u,:) = 0;
        A_cached(u,u) = 1;
        
        u = nz;  % 右下角 (z=Lz, r=0)
        A_cached(u,:) = 0;
        A_cached(u,u) = 1;
        
        % 外部边界 r=Lr
        for i = 1:nz
            u = (nr-1)*nz + i;
            A_cached(u,:) = 0;
            A_cached(u,u) = 1;
        end
        
        % 修正内部点的系数 (r>0)
        for j = 2:nr-1
            for i = 2:nz-1
                u = (j-1)*nz + i;
                r = (j-0.5)*dh;  % 单元中心的r坐标
                
                % 清除现有系数
                A_cached(u,:) = 0;
                
                % 中心点
                A_cached(u,u) = -2/dh^2 - 2/dh^2;
                
                % r方向邻居 (包含1/r项)
                A_cached(u,u-nz) = 1/dh^2 - 1/(2*r*dh);  % r-1
                A_cached(u,u+nz) = 1/dh^2 + 1/(2*r*dh);  % r+1
                
                % z方向邻居
                A_cached(u,u-1) = 1/dh^2;  % z-1
                A_cached(u,u+1) = 1/dh^2;  % z+1
            end
        end
        
        % 添加正则化项以提高数值稳定性
        epsilon = 1e-10;
        A_cached = A_cached + epsilon * speye(size(A_cached));
        
        % 计算LU分解
        [L, U, P, Q] = lu(A_cached);
        last_size = nn;
    end
    
    % 使用迭代求解器代替直接求解
    % 设置迭代参数
    tol = 1e-4;      % 放宽收敛标准
    max_iter = 100;  % 限制迭代次数
    
    % 使用BiCGStab求解器，通常比直接求解更快
    [x_vec, ~] = bicgstab(A_cached, b, tol, max_iter, L, U);
    
    % 如果BiCGStab失败，回退到直接求解
    if isnan(x_vec(1)) || isinf(x_vec(1))
        x_vec = Q * (U \ (L \ (P * b)));
    end
    
    % 重塑回原始网格尺寸
    x = reshape(x_vec, nz, nr);
    
    % 如果使用单精度，将结果转回单精度
    if use_single_precision
        x = single(x);
    end
    
    % 强制边界条件
    x(1,:) = phi0;  % z=0
    x(nz,:) = 0;    % z=Lz
    x(:,nr) = 0;    % r=Lr
    
    % 记录求解时间
    solve_time = toc(tic_local) * 1000;
    if debug_mode && mod(round(w/dh),100) == 0
        fprintf('泊松求解耗时: %.3f ms\n', solve_time);
    end
end


