function A = setup_poisson_matrix(nz, nr, dh, R_matrix, phi0, Lz)
    % 创建泊松方程的系数矩阵 - 优化版本
    % 使用五点差分格式
    % 在圆柱坐标系中，泊松方程为:
    % (1/r)*(d/dr)(r*d(phi)/dr) + d^2(phi)/dz^2 = -rho/epsilon_0
    
    % 声明全局变量
    global use_single_precision
    
    % 计算网格参数
    dr2 = dh^2;
    dz2 = dh^2;
    nn = nz * nr;
    
    % 使用缓存机制提高性能
    persistent A_cached last_nz last_nr last_dh;
    
    % 检查是否可以重用缓存的矩阵
    if ~isempty(A_cached) && last_nz == nz && last_nr == nr && last_dh == dh
        A = A_cached;
        return;  % 直接返回缓存的矩阵，大幅提高性能
    end
    
    % 如果不能重用，构建新矩阵
    tic_matrix = tic;
    
    % 预分配行、列和值数组 - 使用精确的非零元素数量
    nnz_interior = 5 * (nz-2) * (nr-2);  % 内部点的非零元素
    nnz_boundary = 1 * (2*nz + 2*nr - 4);  % 边界点的非零元素
    nnz_total = nnz_interior + nnz_boundary;
    
    rows = zeros(nnz_total, 1);
    cols = zeros(nnz_total, 1);
    vals = zeros(nnz_total, 1);
    
    % 使用向量化操作构建内部点的系数
    count = 0;
    
    % 向量化构建内部点的系数矩阵
    [I, J] = ndgrid(2:nz-1, 2:nr-1);
    indices = sub2ind([nz, nr], I(:), J(:));
    linear_indices = (J(:)-1)*nz + I(:);
    r_values = max(R_matrix(indices), dh/2);  % 避免除以零
    
    % 主对角线系数
    count = count + length(linear_indices);
    rows(1:count) = linear_indices;
    cols(1:count) = linear_indices;
    vals(1:count) = -2/dr2 - 2/dz2;
    
    % z-1 邻居
    start_idx = count + 1;
    end_idx = count + length(linear_indices);
    rows(start_idx:end_idx) = linear_indices;
    cols(start_idx:end_idx) = linear_indices - 1;
    vals(start_idx:end_idx) = 1/dz2;
    count = end_idx;
    
    % z+1 邻居
    start_idx = count + 1;
    end_idx = count + length(linear_indices);
    rows(start_idx:end_idx) = linear_indices;
    cols(start_idx:end_idx) = linear_indices + 1;
    vals(start_idx:end_idx) = 1/dz2;
    count = end_idx;
    
    % r-1 邻居 (包括1/r项)
    start_idx = count + 1;
    end_idx = count + length(linear_indices);
    rows(start_idx:end_idx) = linear_indices;
    cols(start_idx:end_idx) = linear_indices - nz;
    vals(start_idx:end_idx) = 1/dr2 - 1./(2*r_values*dh);
    count = end_idx;
    
    % r+1 邻居 (包括1/r项)
    start_idx = count + 1;
    end_idx = count + length(linear_indices);
    rows(start_idx:end_idx) = linear_indices;
    cols(start_idx:end_idx) = linear_indices + nz;
    vals(start_idx:end_idx) = 1/dr2 + 1./(2*r_values*dh);
    count = end_idx;
    
    % 设置边界条件 - 使用向量化操作
    % 左边界 (z=0)
    left_indices = 1:nz:nn;
    start_idx = count + 1;
    end_idx = count + length(left_indices);
    rows(start_idx:end_idx) = left_indices;
    cols(start_idx:end_idx) = left_indices;
    vals(start_idx:end_idx) = 1.0;
    count = end_idx;
    
    % 右边界 (z=Lz)
    right_indices = nz:nz:nn;
    start_idx = count + 1;
    end_idx = count + length(right_indices);
    rows(start_idx:end_idx) = right_indices;
    cols(start_idx:end_idx) = right_indices;
    vals(start_idx:end_idx) = 1.0;
    count = end_idx;
    
    % 上边界 (r=0)
    top_indices = 1:nz;
    top_indices = setdiff(top_indices, [1, nz]);  % 排除已处理的角点
    
    % 修改r=0轴上的处理，应用L'Hôpital法则
    for i = top_indices
        % 中心点
        rows(count+1) = i;
        cols(count+1) = i;
        vals(count+1) = -4/dh^2 - 2/dh^2;  % 修正中心点系数
        count = count + 1;
        
        % r+1点
        rows(count+1) = i;
        cols(count+1) = i+nz;
        vals(count+1) = 4/dh^2;  % r+1点，系数为4/dh^2
        count = count + 1;
        
        % z-1点
        rows(count+1) = i;
        cols(count+1) = i-1;
        vals(count+1) = 1/dh^2;
        count = count + 1;
        
        % z+1点
        rows(count+1) = i;
        cols(count+1) = i+1;
        vals(count+1) = 1/dh^2;
        count = count + 1;
    end
    
    % 下边界 (r=Lr)
    bottom_indices = (nr-1)*nz+1:nn;
    bottom_indices = setdiff(bottom_indices, [(nr-1)*nz+1, nn]);  % 排除已处理的角点
    start_idx = count + 1;
    end_idx = count + length(bottom_indices);
    rows(start_idx:end_idx) = bottom_indices;
    cols(start_idx:end_idx) = bottom_indices;
    vals(start_idx:end_idx) = 1.0;
    count = end_idx;
    
    % 构建稀疏矩阵
    if use_single_precision
        % 创建稀疏矩阵时必须使用双精度值
        A = sparse(rows(1:count), cols(1:count), double(vals(1:count)), nn, nn);
    else
        A = sparse(rows(1:count), cols(1:count), vals(1:count), nn, nn);
    end
    
    % 缓存矩阵和参数
    A_cached = A;
    last_nz = nz;
    last_nr = nr;
    last_dh = dh;
    
    matrix_time = toc(tic_matrix);
    fprintf('泊松矩阵构建时间: %.3f ms\n', matrix_time*1000);
end