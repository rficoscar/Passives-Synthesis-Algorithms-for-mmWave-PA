function Doherty_Net_Calc_V2()
    % =====================================================================
    % 36GHz mm-Wave Doherty PA Passive Network Synthesis Calculator (V2)
    % Topology: Primary -> 1:n Ideal Xformer -> Shunt Lm -> Series Ll -> Sec
    % =====================================================================

    % 1. 创建主窗口 (UI Figure)
    fig = uifigure('Name', 'mm-Wave Doherty Network Calculator V2', ...
                   'Position',[100, 100, 650, 500]);
               
    % 标题
    uilabel(fig, 'Text', '毫米波 Doherty 功放合成网络参数计算器 (次级Lm-Ll模型)', ...
        'Position',[100, 450, 450, 30], 'FontSize', 16, 'FontWeight', 'bold');

    % ================== 左侧：输入参数区域 ==================
    uilabel(fig, 'Text', '【输入参数 / Inputs】', 'Position',[20, 410, 200, 22], 'FontWeight', 'bold');
    
    uilabel(fig, 'Text', '中心频率 f0 (GHz):', 'Position',[20, 370, 150, 22]);
    ef_f0 = uieditfield(fig, 'numeric', 'Position',[150, 370, 100, 22], 'Value', 36);
    
    uilabel(fig, 'Text', '寄生电容 Cdev (fF):', 'Position',[20, 330, 150, 22]);
    ef_Cdev = uieditfield(fig, 'numeric', 'Position',[150, 330, 100, 22], 'Value', 110);
    
    uilabel(fig, 'Text', '传输线 Z0 (Ohm):', 'Position', [20, 290, 150, 22]);
    ef_Z0 = uieditfield(fig, 'numeric', 'Position',[150, 290, 100, 22], 'Value', 50);
    
    uilabel(fig, 'Text', '传输线长度 Theta (°):', 'Position',[20, 250, 150, 22]);
    ef_theta = uieditfield(fig, 'numeric', 'Position',[150, 250, 100, 22], 'Value', 90);
    
    uilabel(fig, 'Text', '次级并联 Lm (pH):', 'Position',[20, 210, 150, 22]);
    ef_Lm = uieditfield(fig, 'numeric', 'Position',[150, 210, 100, 22], 'Value', 150);
    
    uilabel(fig, 'Text', '次级串联 Ll (pH):', 'Position', [20, 170, 150, 22]);
    ef_Ll = uieditfield(fig, 'numeric', 'Position',[150, 170, 100, 22], 'Value', 50);
    
    uilabel(fig, 'Text', '理想变压器匝比 n:', 'Position',[20, 130, 150, 22]);
    ef_n = uieditfield(fig, 'numeric', 'Position', [150, 130, 100, 22], 'Value', 1.5);

    % ================== 按钮区域 ==================
    btn = uibutton(fig, 'Text', '计算 / Calculate', ...
        'Position',[60, 60, 150, 40], ...
        'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.2 0.6 0.8], 'FontColor', 'w', ...
        'ButtonPushedFcn', @(btn,event) calculateParams());

    % ================== 右侧：输出结果区域 ==================
    uilabel(fig, 'Text', '【输出结果 / Outputs】', 'Position',[320, 410, 200, 22], 'FontWeight', 'bold');
    
    uilabel(fig, 'Text', '1. 抵消 Cdev 所需电感 (pH):', 'Position',[320, 370, 180, 22]);
    out_Lres = uieditfield(fig, 'numeric', 'Position',[500, 370, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    
    uilabel(fig, 'Text', '2. T型 (LCL) 串联 L_T (pH):', 'Position',[320, 330, 180, 22]);
    out_LT = uieditfield(fig, 'numeric', 'Position',[500, 330, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    uilabel(fig, 'Text', '   T型 (LCL) 并联 C_T (fF):', 'Position',[320, 290, 180, 22]);
    out_CT = uieditfield(fig, 'numeric', 'Position',[500, 290, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    
    uilabel(fig, 'Text', '3. Pi型 (CLC) 并联 C_pi (fF):', 'Position',[320, 250, 180, 22]);
    out_Cpi = uieditfield(fig, 'numeric', 'Position',[500, 250, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    uilabel(fig, 'Text', '   Pi型 (CLC) 串联 L_pi (pH):', 'Position',[320, 210, 180, 22]);
    out_Lpi = uieditfield(fig, 'numeric', 'Position', [500, 210, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    
    uilabel(fig, 'Text', '4. 变压器初级电感 Lp (pH):', 'Position',[320, 170, 180, 22]);
    out_Lp = uieditfield(fig, 'numeric', 'Position',[500, 170, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    uilabel(fig, 'Text', '   变压器次级电感 Ls (pH):', 'Position',[320, 130, 180, 22]);
    out_Ls = uieditfield(fig, 'numeric', 'Position',[500, 130, 100, 22], 'Editable', 'off', 'FontColor', 'blue');
    uilabel(fig, 'Text', '   耦合系数 k:', 'Position',[320, 90, 180, 22]);
    out_k = uieditfield(fig, 'numeric', 'Position',[500, 90, 100, 22], 'Editable', 'off', 'FontColor', 'blue');

    % 初始化计算一次
    calculateParams();

    % ================== 回调计算函数 ==================
    function calculateParams()
        % 提取输入参数
        f0_GHz  = ef_f0.Value;
        Cdev_fF = ef_Cdev.Value;
        Z0      = ef_Z0.Value;
        theta   = ef_theta.Value;
        Lm_pH   = ef_Lm.Value;
        Ll_pH   = ef_Ll.Value;
        n       = ef_n.Value;
        
        % 物理常数与单位换算
        w0 = 2 * pi * f0_GHz * 1e9;
        Cdev = Cdev_fF * 1e-15;
        theta_rad = theta * pi / 180;
        
        % 1. 计算谐振电感 L_res
        L_res = 1 / (w0^2 * Cdev);
        out_Lres.Value = L_res * 1e12; % 转 pH
        
        % 2. 传输线 T 型等效计算
        L_T = Z0 * tan(theta_rad / 2) / w0;
        C_T = sin(theta_rad) / (w0 * Z0);
        out_LT.Value = L_T * 1e12; % 转 pH
        out_CT.Value = C_T * 1e15; % 转 fF
        
        % 3. 传输线 Pi 型等效计算
        C_pi = tan(theta_rad / 2) / (w0 * Z0);
        L_pi = Z0 * sin(theta_rad) / w0;
        out_Cpi.Value = C_pi * 1e15; % 转 fF
        out_Lpi.Value = L_pi * 1e12; % 转 pH
        
        % 4. 变压器参数换算 (基于 次级Lm, Ll 模型)
        % 【关键修改点：新的拓扑公式】
        Lp = Lm_pH / (n^2);
        Ls = Lm_pH + Ll_pH;
        k  = sqrt(Lm_pH / (Lm_pH + Ll_pH));
        
        out_Lp.Value = Lp;
        out_Ls.Value = Ls;
        out_k.Value  = k;
    end
end