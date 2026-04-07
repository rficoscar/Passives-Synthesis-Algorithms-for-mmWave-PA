function dual_band_optimization_jssc_style
    clc; clear; close all;

    %% 1. 物理参数与目标设定
    f1 = 40e9;              % 频点 1: 40 GHz
    alpha = 25/40;          
    f2 = f1 * alpha;        % 频点 2: 25 GHz
    
    R_target = 50;          % Port 1 源阻抗 (Rp)
    RL = 50;                % Port 2 负载阻抗 (Rs)
    
    % --- 搜索空间限制 ---
    Cp_lim = [100, 100];      % 100fF
    Cs_lim = [100, 100];       % 100fF
    L_lim  = [100, 1000];   % pH
    k_lim  = [0.2, 0.9];  % k

    % 变量: x = [Lp, Ls, k, Cp, Cs]
    lb = [L_lim(1), L_lim(1), k_lim(1), Cp_lim(1), Cs_lim(1)];
    ub = [L_lim(2), L_lim(2), k_lim(2), Cp_lim(2), Cs_lim(2)];
    x0 = (lb + ub) / 2;

    %% 2. 优化设置
    history.x = []; history.fval = []; history.iteration = [];
    
    options = optimoptions('fmincon', ...
        'Display', 'iter-detailed', ... 
        'Algorithm', 'sqp', ...         
        'StepTolerance', 1e-10, ...
        'MaxFunctionEvaluations', 2000, ...
        'OutputFcn', @outfun);          

    %% 3. 执行优化
    fprintf('开始优化...\n');
    tic;
    [x_opt, fval] = fmincon(@objective_func, x0, [],[],[],[], lb, ub, [], options);
    solve_time = toc;

    %% 4. 结果提取
    Lp_final = x_opt(1) * 1e-12;
    Ls_final = x_opt(2) * 1e-12;
    k_final  = x_opt(3);
    Cp_final = x_opt(4) * 1e-15;
    Cs_final = x_opt(5) * 1e-15;

    fprintf('\n========== 优化结果 ==========\n');
    fprintf('耗时: %.2f 秒\n', solve_time);
    fprintf('Lp = %.2f pH, Ls = %.2f pH, k = %.4f\n', x_opt(1), x_opt(2), x_opt(3));
    fprintf('Cp = %.2f fF, Cs = %.2f fF\n', x_opt(4), x_opt(5));
    
    %% 5. 生成 JSSC 风格验证图 (新窗口)
    verify_result_in_new_window(Lp_final, Ls_final, k_final, Cp_final, Cs_final, R_target, RL, f1, f2);


    %% --- 内部函数定义 ---

    % 1. 目标函数
    function err = objective_func(x)
        Lp = x(1) * 1e-12; Ls = x(2) * 1e-12; k_v = x(3);
        Cp = x(4) * 1e-15; Cs = x(5) * 1e-15;
        Lm = k_v * sqrt(Lp * Ls);
        
        Z1 = calc_Zin_point(f1, Lp, Ls, Lm, Cp, Cs, k_v, RL);
        Z2 = calc_Zin_point(f2, Lp, Ls, Lm, Cp, Cs, k_v, RL);
        
        w_real = 1.0; w_imag = 1.5; 
        err1 = w_real*((real(Z1) - R_target)/R_target)^2 + w_imag*(imag(Z1)/R_target)^2;
        err2 = w_real*((real(Z2) - R_target)/R_target)^2 + w_imag*(imag(Z2)/R_target)^2;
        err = err1 + err2;
    end

    % 2. 单点 Zin 计算
    function Zin = calc_Zin_point(f, Lp, Ls, Lm, Cp, Cs, k_val, RL_val)
        w = 2 * pi * f; j = 1i;
        Y11 = j*w*Cp + 1/(j*w*Lp*(1-k_val^2));
        Y22 = j*w*Cs + 1/(j*w*Ls*(1-k_val^2));
        Y12 = -k_val^2 / (j*w*Lm*(1-k_val^2));
        Y21 = Y12;
        Yin = Y11 - (Y12*Y21)/(Y22 + 1/RL_val);
        Zin = 1/Yin;
    end

    % 3. 单点 S参数计算 (Port1=Rp, Port2=Rs)
    function [S11, S21] = calc_S_params(f, Lp, Ls, k_val, Cp, Cs, Rp, Rs)
        w = 2 * pi * f; j = 1i;
        Lm = k_val * sqrt(Lp * Ls);
        
        % Y矩阵
        Y11 = j*w*Cp + 1/(j*w*Lp*(1-k_val^2));
        Y22 = j*w*Cs + 1/(j*w*Ls*(1-k_val^2));
        Y12 = -k_val^2 / (j*w*Lm*(1-k_val^2));
        Y21 = Y12;
        
        % S参数转换 (Power Wave formulation for arbitrary impedances)
        gp = 1/Rp; % Port 1 conductance
        gs = 1/Rs; % Port 2 conductance
        
        % 特征分母
        Delta = (Y11 + gp)*(Y22 + gs) - Y12*Y21;
        
        % S21 (Transducer Gain相关)
        S21 = -2 * Y21 * sqrt(gp * gs) / Delta;
        
        % S11
        % 先算输入导纳 Yin_loaded
        Yin_loaded = Y11 - (Y12*Y21)/(Y22 + gs);
        % 反射系数
        S11 = (gp - Yin_loaded) / (gp + Yin_loaded); % 注意导纳形式的反射系数公式差异
        % 或者用阻抗形式更直观: Zin = 1/Yin_loaded; S11 = (Zin - Rp)/(Zin + Rp);
        Zin_val = 1/Yin_loaded;
        S11 = (Zin_val - Rp) / (Zin_val + Rp);
    end

    % 4. 优化过程可视化
    function stop = outfun(x, optimValues, state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
                history.fval = [history.fval; optimValues.fval];
                history.iteration = [history.iteration; optimValues.iteration];
                hFig = findobj('Type','figure','Tag','LivePlot');
                if isempty(hFig)
                    hFig = figure('Name','Optimization Progress', 'Tag', 'LivePlot', ...
                                  'Color', 'w', 'Position', [100, 200, 1000, 300]); 
                end
                set(0, 'CurrentFigure', hFig); 
                subplot(1,3,1); plot(history.iteration, history.x(:,1), 'b-', history.iteration, history.x(:,2), 'r--'); xlabel('Iteration Number');
                title('(a) Inductance (pH)'); legend('Lp','Ls'); grid on; xlim([0, max(1,optimValues.iteration)]);
                subplot(1,3,2); plot(history.iteration, history.x(:,3), 'k-'); xlabel('Iteration Number');
                title('(b) Coupling k'); grid on; xlim([0, max(1,optimValues.iteration)]); ylim([0,1]);
                subplot(1,3,3); semilogy(history.iteration, history.fval, 'm-'); xlabel('Iteration Number');
                title('(c) Error'); grid on; xlim([0, max(1,optimValues.iteration)]);
                drawnow limitrate;
        end
    end

    % 5. 最终验证绘图 (JSSC 风格)
    function verify_result_in_new_window(Lp, Ls, k, Cp, Cs, Rp, Rs, f1, f2)
        % 扫描范围
        freq = linspace(10e9, 60e9, 501);
        Z_real = zeros(size(freq));
        Z_imag = zeros(size(freq));
        S11_dB = zeros(size(freq));
        S21_dB = zeros(size(freq));
        Lm = k * sqrt(Lp * Ls);
        
        for i = 1:length(freq)
            % 计算 Zin
            z = calc_Zin_point(freq(i), Lp, Ls, Lm, Cp, Cs, k, Rs);
            Z_real(i) = real(z);
            Z_imag(i) = imag(z);
            
            % 计算 S参数 (Port1=Rp, Port2=Rs)
            [s11, s21] = calc_S_params(freq(i), Lp, Ls, k, Cp, Cs, Rp, Rs);
            S11_dB(i) = 20*log10(abs(s11));
            S21_dB(i) = 20*log10(abs(s21));
        end
        
        % 创建 JSSC 风格窗口
        hJSSC = figure('Name', 'Final Verification (JSSC Style)', 'Color', 'w', 'Position', [150, 150, 900, 400]);
        
        % --- 子图 1: Zin ---
        subplot(1, 2, 1);
        plot(freq/1e9, Z_real, 'r-', 'LineWidth', 1.5); hold on; % 实部用黑色实线
        plot(freq/1e9, Z_imag, 'b-', 'LineWidth', 1.5);        % 虚部用蓝色虚线
        
        % 标记频点
        xline(f1/1e9, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
        xline(f2/1e9, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
        yline(0, 'k-', 'LineWidth', 0.5); % 0线
		yline(Rp, 'k--', 'LineWidth', 0.5);
        
        % 样式修饰
        title('(a) Input Impedance');
        xlabel('Frequency (GHz)');
        ylabel('Impedance (\Omega)');
        legend('Re(Z_{in})', 'Im(Z_{in})', 'Location', 'best');
        xlim([min(freq)/1e9, max(freq)/1e9]);
        % ylim([-100, 150]); % 根据需要调整
        
        % JSSC 字体与边框设置
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1, 'Box', 'on');
        grid on;
        
        % --- 子图 2: S-Parameters ---
        subplot(1, 2, 2);
        plot(freq/1e9, S21_dB, 'r-', 'LineWidth', 1.5); hold on; % S21 用黑色实线
        plot(freq/1e9, S11_dB, 'b-', 'LineWidth', 1.5);        % S11 用蓝色虚线
        
        % 标记频点与参考线
        xline(f1/1e9, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
        xline(f2/1e9, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);
        yline(-10, 'r:', 'LineWidth', 1); % -10dB 参考线
        
        % 样式修饰
        title('(b) S-Parameters');
        xlabel('Frequency (GHz)');
        ylabel('Magnitude (dB)');
        legend('S_{21}', 'S_{11}', 'Location', 'best');
        xlim([min(freq)/1e9, max(freq)/1e9]);
        ylim([-40, 5]); % 典型 S 参数显示范围
        
        % JSSC 字体与边框设置
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1, 'Box', 'on');
        grid on;
        
        % 标注 Port 阻抗信息
        % sgtitle(['Transformer Network Analysis (Port1 Z_S=', num2str(Rp), '\Omega, Port2 Z_L=', num2str(Rs), '\Omega)'], ...
        %         'FontName', 'Times New Roman', 'FontSize', 14);
    end

end