function gain_slope_compensation
    clc; clear; close all;

    %% 1. 系统参数设置
    % 频率设置
    f1 = 40e9;              % 高频点 40 GHz
    f2 = 25e9;              % 低频点 25 GHz
    
    % 端口阻抗
    R_source = 50;          % Port 1 (源) 阻抗 (对应你之前的 Rtarget)
    R_load   = 50;          % Port 2 (负载) 阻抗 (对应 RL)

    % 晶体管增益参数 (线性模型)
    % 假设在 25GHz 到 40GHz 之间下降 3.5dB (取3和4的平均值)
    Gain_drop = 3.5;        % dB
    G_Tx_ref  = 10;         % 晶体管在 f2 (25GHz) 处的参考增益 (dB) - 绝对值不影响斜率优化，只影响绘图高度
    
    % 定义晶体管增益函数句柄 (单位: GHz -> dB)
    % 斜率 k = -Drop / (f1 - f2)
    slope_Tx = -Gain_drop / ((f1 - f2)/1e9); 
    get_Tx_Gain = @(f_hz) G_Tx_ref + slope_Tx * (f_hz/1e9 - f2/1e9);

    %% 2. 搜索范围 (pH, fF)
    % 物理单位: pH, fF. 
    Cp_lim = [70, 70];      % fF
    Cs_lim = [70, 70];      % fF
    L_lim  = [50, 800];     % pH (稍微放宽下限，以获得更好的高频特性)
    k_lim  = [0.2, 0.9];  % k

    % 变量归一化: x = [Lp, Ls, k, Cp, Cs]
    lb = [L_lim(1), L_lim(1), k_lim(1), Cp_lim(1), Cs_lim(1)];
    ub = [L_lim(2), L_lim(2), k_lim(2), Cp_lim(2), Cs_lim(2)];
    x0 = (lb + ub) / 2;

    %% 3. 优化设置 (加入可视化)
    history.fval = [];
    history.iter = [];
    
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ...
        'Algorithm', 'sqp', ...
        'StepTolerance', 1e-10, ...
        'MaxFunctionEvaluations', 3000, ...
        'OutputFcn', @outfun);

    fprintf('开始增益斜率补偿优化 (Gain Slope Compensation)...\n');
    fprintf('目标: 补偿晶体管 %.1f dB 的增益滚降\n', Gain_drop);
    
    [x_opt, fval] = fmincon(@cost_function, x0, [],[],[],[], lb, ub, [], options);

    %% 4. 结果分析与绘图
    % 还原参数
    Lp = x_opt(1)*1e-12; Ls = x_opt(2)*1e-12; k = x_opt(3);
    Cp = x_opt(4)*1e-15; Cs = x_opt(5)*1e-15;
    
    fprintf('\n====== 优化完成 ======\n');
    fprintf('Lp = %.2f pH\n', x_opt(1));
    fprintf('Ls = %.2f pH\n', x_opt(2));
    fprintf('k  = %.4f\n', x_opt(3));
    fprintf('Cp = %.2f fF\n', x_opt(4));
    fprintf('Cs = %.2f fF\n', x_opt(5));
    
    % 最终全频段扫描
    f_sweep = linspace(10e9, 60e9, 501);
    
    % 计算数组
    S21_net_dB = zeros(size(f_sweep));
    S11_net_dB = zeros(size(f_sweep));
    G_Tx_dB    = zeros(size(f_sweep));
    
    for i = 1:length(f_sweep)
        [s11, s21] = calc_S_params(f_sweep(i), Lp, Ls, k, Cp, Cs, R_source, R_load);
        S21_net_dB(i) = 20*log10(abs(s21));
        S11_net_dB(i) = 20*log10(abs(s11));
        G_Tx_dB(i)    = get_Tx_Gain(f_sweep(i));
    end
    
    G_Total_dB = S21_net_dB + G_Tx_dB;
    
    % --- 绘制结果对比图 ---
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    
    % 1. 增益分析
    subplot(2, 1, 1);
    plot(f_sweep/1e9, G_Tx_dB, 'k--', 'LineWidth', 1.5); hold on;
    plot(f_sweep/1e9, S21_net_dB, 'b-', 'LineWidth', 1.5);
    plot(f_sweep/1e9, G_Total_dB, 'r-', 'LineWidth', 2.5);
    
    xline(f1/1e9, 'k:', 'LineWidth', 1);
    xline(f2/1e9, 'k:', 'LineWidth', 1);
    
    title('Gain Slope Compensation');
    xlabel('Frequency (GHz)');
    ylabel('Gain (dB)');
    legend('Transistor Gain (Model)', 'Matching Network S21', 'Total Gain (Flat Target)', 'Location', 'best');
    grid on;
    
    % 标注关键点
    idx_f1 = find(f_sweep >= f1, 1);
    idx_f2 = find(f_sweep >= f2, 1);
    if ~isempty(idx_f2), text(f2/1e9, G_Total_dB(idx_f2)+1, sprintf('%.2fdB', G_Total_dB(idx_f2)), 'Color', 'r'); end
    if ~isempty(idx_f1), text(f1/1e9, G_Total_dB(idx_f1)+1, sprintf('%.2fdB', G_Total_dB(idx_f1)), 'Color', 'r'); end

    % 2. S11 (虽然不是优化目标，但需要看一眼)
    subplot(2, 1, 2);
    plot(f_sweep/1e9, S11_net_dB, 'm-', 'LineWidth', 1.5);
    yline(-10, 'k--');
    title('Input Return Loss (S11)');
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    ylim([-30, 0]);


    %% --- 内部函数 ---

    % 代价函数
    function err = cost_function(x)
        Lp_v = x(1)*1e-12; Ls_v = x(2)*1e-12; k_v = x(3);
        Cp_v = x(4)*1e-15; Cs_v = x(5)*1e-15;
        
        % 计算 f1 和 f2 处的 S21
        [~, s21_f1] = calc_S_params(f1, Lp_v, Ls_v, k_v, Cp_v, Cs_v, R_source, R_load);
        [~, s21_f2] = calc_S_params(f2, Lp_v, Ls_v, k_v, Cp_v, Cs_v, R_source, R_load);
        
        S21_dB_f1 = 20*log10(abs(s21_f1));
        S21_dB_f2 = 20*log10(abs(s21_f2));
        
        % 晶体管增益
        G_Tx_f1 = get_Tx_Gain(f1);
        G_Tx_f2 = get_Tx_Gain(f2);
        
        % 总增益
        Total_f1 = G_Tx_f1 + S21_dB_f1;
        Total_f2 = G_Tx_f2 + S21_dB_f2;
        
        % 目标1: 平坦度 (差值平方)
        flatness_error = (Total_f1 - Total_f2)^2;
        
        % 目标2: 插入损耗惩罚 (避免优化器找出 S21 = -100dB 的平坦解)
        % 我们希望 S21 尽可能大（接近 0dB）
        loss_penalty = 0.1 * (abs(S21_dB_f1) + abs(S21_dB_f2)); 
        
        err = flatness_error + loss_penalty;
    end

    % 计算 S参数核心函数
    function [S11, S21] = calc_S_params(f, Lp, Ls, k, Cp, Cs, Rs, Rl)
        w = 2 * pi * f;
        j = 1i;
        Lm = k * sqrt(Lp * Ls);
        
        % Y参数计算 (Cps = 0)
        Y11 = j*w*Cp + 1/(j*w*Lp*(1-k^2));
        Y22 = j*w*Cs + 1/(j*w*Ls*(1-k^2));
        Y12 = -k^2 / (j*w*Lm*(1-k^2));
        Y21 = Y12;
        
        % Y -> S 参数转换 (针对任意实数参考阻抗 Rs, Rl)
        % 定义参考导纳
        Gs = 1/Rs; 
        Gl = 1/Rl;
        
        % 分母 Delta
        Delta = (Y11 + Gs)*(Y22 + Gl) - Y12*Y21;
        
        % S21 公式 (Transducer Power Gain related)
        % S21 = -2 * Y21 * sqrt(Gs * Gl) / Delta
        S21 = -2 * Y21 * sqrt(Gs * Gl) / Delta;
        
        % S11 公式
        % Yin = Y11 - Y12*Y21 / (Y22 + Gl)
        Yin = Y11 - (Y12*Y21) / (Y22 + Gl);
        % Zin = 1/Yin
        % Gamma = (Zin - Rs) / (Zin + Rs)
        % 或者直接用 Y 参数计算 S11
        S11 = ((1/Yin) - Rs) / ((1/Yin) + Rs);
    end

    % 可视化回调
    function stop = outfun(x, optimValues, state)
        stop = false;
        if strcmp(state, 'iter')
            history.fval = [history.fval; optimValues.fval];
            history.iter = [history.iter; optimValues.iteration];
            
            set(0, 'CurrentFigure', findobj('Type','figure','Tag','LivePlot'));
            if isempty(findobj('Type','figure','Tag','LivePlot'))
                figure('Name','Opt Process', 'Tag', 'LivePlot', 'Color', 'w', 'Position', [50, 50, 800, 400]);
            end
            
            subplot(1, 2, 1);
            plot(history.iter, history.fval, 'b.-');
            title('Cost Convergence'); xlabel('Iter'); grid on; set(gca, 'YScale', 'log');
            
            subplot(1, 2, 2);
            bar([x(1), x(2); x(4), x(5)]);
            title('Current L(pH) & C(fF)'); legend('Pri', 'Sec');
            xticklabels({'L', 'C'});
            drawnow;
        end
    end
end