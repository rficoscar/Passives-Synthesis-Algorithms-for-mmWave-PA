function solve_impedance_final()
    clc; clear; close all;

    %% 1. 物理参数与配置
    % --- 频率设置 ---
    f0 = 40e9; w0 = 2 * pi * f0;
    f1 = 25e9; w1 = 2 * pi * f1;
    alpha = w1/w0;
    
    % --- 基础器件参数 ---
    GL = 0.02;              
    Cdev_base = 260e-15;    % 260 fF
    Gopt_base = 0.04;       
    
    % --- 扫描范围 ---
    cl_range_val = [0, 200e-15];    % CL 范围
    scale_range_val = [0.5, 1.5];   % Scale 范围
    theta_range_val = [10, 60];     % Theta 范围
    
    % --- 误差权重 (用户要求: 提高 Zo 权重) ---
    W_zo = 5.0;  % Zo 误差的权重
    W_ze = 1.0;  % Ze 误差的权重
    
    % --- 迭代策略 ---
    target_solutions = 5;    % 至少找几个解
    time_limit = 60;         % 超时时间(秒)
    tolerance = 1e-1;        % 误差容限
    
    % 初始网格密度
    n_scale = 1000; 
    n_cl = 100;
    max_scale_pts = 1000;
    max_cl_pts = 100;

    %% 2. 迭代求解循环
    t_start = tic;
    iter_info = []; % 存储迭代日志用于MD输出
    
    % 最终结果容器
    final_X_CL = []; final_Y_Scale = [];
    final_Error_Map = [];
    final_C1_Map = []; final_C2_Map = [];
    final_solutions = []; 
    
    fprintf('开始计算 (公式已修正, Zo权重=%.1f)...\n', W_zo);
    
    while true
        % 2.1 生成网格
        scale_vec = linspace(scale_range_val(1), scale_range_val(2), n_scale);
        cl_vec = linspace(cl_range_val(1), cl_range_val(2), n_cl);
        
        % 预分配
        map_err = nan(n_scale, n_cl);
        map_theta = zeros(n_scale, n_cl);
        map_c1 = false(n_scale, n_cl);
        map_c2 = false(n_scale, n_cl);
        
        % 超时检查
        if toc(t_start) > time_limit, break; end
        
        % 2.2 网格扫描
        for i = 1:n_scale
            s = scale_vec(i);
            Cdev = Cdev_base * s;
            Gopt = Gopt_base * s;
            
            for j = 1:n_cl
                CL = cl_vec(j);
                
                % (A) 优化 Theta 最小化加权误差
                obj = @(t) calc_weighted_error(t, w0, alpha, Cdev, CL, Gopt, GL, W_zo, W_ze);
                [best_rad, min_err] = fminbnd(obj, theta_range_val(1)*pi/180, theta_range_val(2)*pi/180);
                
                map_err(i,j) = min_err;
                map_theta(i,j) = best_rad * 180/pi;
                
                % (B) 检查物理约束 (使用修正后的公式)
                % Set 1
                [~,~,yp1,ym1] = calc_imp_corrected(w0, best_rad, Cdev, CL, Gopt, GL);
                map_c1(i,j) = (yp1 > ym1) && (ym1 > 0);
                
                % Set 2
                [~,~,yp2,ym2] = calc_imp_corrected(alpha*w0, alpha*best_rad, Cdev, CL, Gopt, GL);
                map_c2(i,j) = (yp2 > ym2) && (ym2 > 0);
            end
        end
        
        % 2.3 统计解
        valid_mask = (map_err < tolerance) & map_c1 & map_c2;
        num_found = sum(valid_mask(:));
        
        % 记录日志
        log_str = sprintf('Iter: ScalePts=%d, CLPts=%d, Time=%.1fs, Found=%d', ...
            n_scale, n_cl, toc(t_start), num_found);
        fprintf('%s\n', log_str);
        iter_info{end+1} = log_str; %#ok<AGROW>
        
        % 保存当前数据用于绘图
        [final_X_CL, final_Y_Scale] = meshgrid(cl_vec, scale_vec);
        final_Error_Map = map_err;
        final_C1_Map = map_c1;
        final_C2_Map = map_c2;
        final_Theta_Map = map_theta;
        
        % 2.4 判断退出
        if num_found >= target_solutions || toc(t_start) > time_limit
            break;
        end
        
        % 2.5 增加精度策略 (优先Scale)
        if n_scale < max_scale_pts
            n_scale = n_scale * 2;
        elseif n_cl < max_cl_pts
            n_cl = n_cl * 2;
        else
            break; 
        end
    end
    
    %% 3. 提取所有解的详细信息
    [r_idx, c_idx] = find((final_Error_Map < tolerance) & final_C1_Map & final_C2_Map);
    
    fprintf('正在整理 %d 个解的详细数据...\n', length(r_idx));
    sol_struct = [];
    
    for k = 1:length(r_idx)
        r = r_idx(k); c = c_idx(k);
        
        s = final_Y_Scale(r,c);
        cl_val = final_X_CL(r,c);
        theta_deg = final_Theta_Map(r,c);
        theta_rad = theta_deg * pi/180;
        
        cdev_val = Cdev_base * s;
        gopt_val = Gopt_base * s;
        
        % 计算两组频率下的详细值
        [z1o, z1e, yp1, ym1] = calc_imp_corrected(w0, theta_rad, cdev_val, cl_val, gopt_val, GL);
        [z2o, z2e, yp2, ym2] = calc_imp_corrected(alpha*w0, alpha*theta_rad, cdev_val, cl_val, gopt_val, GL);
        
        % 存入结构体
        sol_struct(k).id = k;
        sol_struct(k).scale = s;
        sol_struct(k).cl = cl_val;
        sol_struct(k).cdev = cdev_val;
        sol_struct(k).gopt = gopt_val;
        sol_struct(k).theta = theta_deg;
        sol_struct(k).error = final_Error_Map(r,c);
        
        sol_struct(k).z1o = z1o; sol_struct(k).z1e = z1e;
        sol_struct(k).yp1 = yp1; sol_struct(k).ym1 = ym1;
        
        sol_struct(k).z2o = z2o; sol_struct(k).z2e = z2e;
        sol_struct(k).yp2 = yp2; sol_struct(k).ym2 = ym2;
    end

    %% 4. 绘图与导出 SVG
    
    % --- Figure 1: 误差分布与最终解 (2 subplots) ---
    f1_handle = figure('Position', [100, 100, 1200, 500], 'Color', 'w');
    
    % Subplot 1: 误差热力图
    subplot(1, 2, 1);
    surf(final_X_CL*1e15, final_Y_Scale, log10(final_Error_Map+1e-12));
    view(2); shading interp; colorbar;
    title('Log10(Weighted Error)');
    xlabel('CL (fF)'); ylabel('Scale');
    pbaspect([1 0.7 1]); axis tight;
    
    % Subplot 2: 最终解三维分布
    subplot(1, 2, 2);
    if ~isempty(sol_struct)
        scatter3([sol_struct.cl]*1e15, [sol_struct.scale], [sol_struct.theta], ...
            40, [sol_struct.error], 'filled');
        colorbar;
        title(sprintf('Valid Solutions (Total: %d)', length(sol_struct)));
        xlabel('CL (fF)'); ylabel('Scale'); zlabel('Theta (deg)');
        grid on; view(3);
    else
        text(0.5,0.5,'No Solution','HorizontalAlignment','center');
    end
    pbaspect([1 0.7 1]);
    
    % 导出 Figure 1
    saveas(f1_handle, 'Result_Analysis.svg');
    fprintf('已导出: Result_Analysis.svg\n');
    
    % --- Figure 2: Yp Ym 约束检查 (1行3列) ---
    f2_handle = figure('Position', [100, 100, 1500, 400], 'Color', 'w');
    
    % Subplot 1: Set 1 Constraints
    subplot(1, 3, 1);
    % 使用 imagesc 绘图需要注意坐标轴方向
    surf(final_X_CL*1e15, final_Y_Scale, double(final_C1_Map));
    view(2); shading flat; colormap(gca, [0.2 0.2 0.8; 1 1 0]); % 蓝=False, 黄=True
    title('Constraint: Yp1 > Ym1 > 0');
    xlabel('CL (fF)'); ylabel('Scale');
    axis tight; pbaspect([1 0.8 1]);
    
    % Subplot 2: Set 2 Constraints
    subplot(1, 3, 2);
    surf(final_X_CL*1e15, final_Y_Scale, double(final_C2_Map));
    view(2); shading flat; colormap(gca, [0.2 0.2 0.8; 1 1 0]);
    title('Constraint: Yp2 > Ym2 > 0');
    xlabel('CL (fF)'); ylabel('Scale');
    axis tight; pbaspect([1 0.8 1]);
    
    % Subplot 3: Intersection
    subplot(1, 3, 3);
    surf(final_X_CL*1e15, final_Y_Scale, double(final_C1_Map & final_C2_Map));
    view(2); shading flat; colormap(gca, [0.2 0.2 0.8; 1 1 0]);
    title('Intersection (Both Valid)');
    xlabel('CL (fF)'); ylabel('Scale');
    axis tight; pbaspect([1 0.8 1]);
    
    % 导出 Figure 2
    saveas(f2_handle, 'Constraints_Check.svg');
    fprintf('已导出: Constraints_Check.svg\n');

    %% 5. 导出 Markdown 报告
    write_markdown('Impedance_Report.md', sol_struct, iter_info, ...
        f0, f1, alpha, GL, Cdev_base, Gopt_base, W_zo, W_ze);

end

%% --- 辅助计算函数 (公式修正版) ---

function e = calc_weighted_error(t, w, a, Cdev, CL, Gopt, GL, w_zo, w_ze)
    [z1o, z1e, ~, ~] = calc_imp_corrected(w, t, Cdev, CL, Gopt, GL);
    [z2o, z2e, ~, ~] = calc_imp_corrected(a*w, a*t, Cdev, CL, Gopt, GL);
    
    % 加权误差计算
    err_zo = abs(z1o - z2o);
    err_ze = abs(z1e - z2e);
    
    e = w_zo * err_zo + w_ze * err_ze;
end

function [z0o, z0e, yp, ym] = calc_imp_corrected(w, t, Cdev, CL, Gopt, GL)
    % 1. 计算 Yp
    k = Gopt / GL; % 用于 Yp 公式中的比率
    num = w * Cdev - k * w * CL;
    den = cot(t) - k * cot(2*t);
    
    if abs(den) < 1e-9, den = 1e-9; end
    yp = num / den;
    
    % 2. 计算 Ym (公式修正处)
    % 图片公式: Ym = sin(t) * sqrt( (2*Gopt/GL) * [GL^2 + (Yp*cot2t - wCL)^2] )
    % 旧代码错误: term = (2*k/GL) ... 其中 k=Gopt/GL -> 2*Gopt/GL^2 (错)
    % 修正代码:
    factor = 2 * Gopt / GL; 
    
    inner_term = GL^2 + (yp * cot(2*t) - w * CL)^2;
    val_in_sqrt = factor * inner_term;
    
    % 安全开根号 (防止负数报错, 物理约束会在外部检查)
    ym = sin(t) * sqrt(abs(val_in_sqrt));
    
    % 3. 计算 Z0
    z0o = 1 / (yp + ym);
    z0e = 1 / (yp - ym);
end

%% --- Markdown 写入函数 ---

function write_markdown(fname, sols, logs, f0, f1, alpha, GL, Cdev_b, Gopt_b, wzo, wze)
    fid = fopen(fname, 'w', 'n', 'UTF-8');
    if fid == -1, return; end
    
    fprintf('正在写入文件 %s ...\n', fname);
    
    fprintf(fid, '# Impedance Simulation Report\n\n');
    fprintf(fid, '**Date**: %s\n\n', datestr(now));
    
    fprintf(fid, '## 1. Simulation Settings\n');
    fprintf(fid, '- **Frequencies**: %.2f GHz / %.2f GHz (Alpha=%.4f)\n', f0/1e9, f1/1e9, alpha);
    fprintf(fid, '- **Base Params**: Cdev=%.0f fF, Gopt=%.3f S, GL=%.3f S\n', Cdev_b*1e15, Gopt_b, GL);
    fprintf(fid, '- **Error Weights**: Zo_Weight = **%.1f**, Ze_Weight = %.1f\n', wzo, wze);
    fprintf(fid, '- **Formula Fix**: Ym calculation uses `2*Gopt/GL` pre-factor.\n\n');
    
    fprintf(fid, '## 2. Iteration Log\n');
    fprintf(fid, 'Iterative grid refinement process:\n\n');
    fprintf(fid, '```text\n');
    for i = 1:length(logs)
        fprintf(fid, '%s\n', logs{i});
    end
    fprintf(fid, '```\n\n');
    
    fprintf(fid, '## 3. Solution Data Table\n');
    if isempty(sols)
        fprintf(fid, 'No valid solutions found.\n');
    else
        fprintf(fid, 'Detailed data for all valid solutions:\n\n');
        
        % 表头
        fprintf(fid, '| ID | Scale | CL(fF) | Cdev(fF) | Gopt(S) | Theta(deg) | Error | Z1o | Z2o | Yp1(S) | Ym1(S) | Yp2(S) | Ym2(S) |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
        
        for k = 1:length(sols)
            s = sols(k);
            % 这里列出关键信息，由于列数太多，Z1e/Z2e 暂时省略，或者合并显示
            fprintf(fid, '| %d | %.4f | %.2f | %.2f | %.4f | %.2f | %.2e | %.1f | %.1f | %.2e | %.2e | %.2e | %.2e |\n', ...
                s.id, s.scale, s.cl*1e15, s.cdev*1e15, s.gopt, s.theta, s.error, ...
                s.z1o, s.z2o, ...
                s.yp1, s.ym1, s.yp2, s.ym2);
        end
    end
    
    fclose(fid);
    fprintf('MD 报告已生成.\n');
end