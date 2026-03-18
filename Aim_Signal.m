%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 日期：2022/3/6（修正版：2025）
% 函数名称：生成矢量水听器接收目标信号
% 功能：双目标 + 多线谱 + 有色背景 + 白噪声干扰
% 【重要更新】SNR 精确控制：信号与环境噪声同步通过FIR滤波器
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data] = Aim_Signal(N, Num_N, t, frequency, theta, SNR, p_sx, Fir)

    Data = zeros(4, N * Num_N);
    deread = pi / 180;
    
    % 预计算滤波器（-6dB 噪声模型）：用于生成目标自身有色背景
    SOS = [1, 1, 0, 1, -0.947134927964088, 0];
    G = 0.026432536017956;
    [b_color, a_color] = sos2tf(SOS, G);

    for iiii = 1:Num_N
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. 生成每个目标的“辐射噪声”（含线谱 + 有色背景）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % --- 目标1 ---
        p1_line = p_sx(1)*cos(2*pi*frequency(1)*t) + ...
                  p_sx(2)*cos(2*pi*frequency(2)*t) + ...
                  0.25*cos(2*pi*90*t)  + 0.21*cos(2*pi*130*t) + ...
                  0.28*cos(2*pi*170*t)   + 0.24*cos(2*pi*120*t);
                  
        background1 = filter(b_color, a_color, randn(size(t)));
        p1 = p1_line + background1;  % 背景作为目标自身辐射噪声
        
        vx1 = p1 * cos(theta(1) * deread);
        vy1 = p1 * sin(theta(1) * deread);
        
        % --- 目标2 ---
        p2_line = p_sx(3)*cos(2*pi*frequency(3)*t) + ...
                  p_sx(4)*cos(2*pi*frequency(4)*t) + ...
                  0.28*cos(2*pi*110*t)  + 0.26*cos(2*pi*140*t) + ...
                  0.29*cos(2*pi*150*t)  + 0.27*cos(2*pi*160*t);
                  
        background2 = filter(b_color, a_color, randn(size(t)));
        p2 = p2_line + background2;
        
        vx2 = p2 * cos(theta(2) * deread);
        vy2 = p2 * sin(theta(2) * deread);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. 合成总干净信号（无外部噪声）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_clean = p1 + p2;
        vx_clean = vx1 + vx2;
        vy_clean = vy1 + vy2;
        vz_clean = zeros(size(p_clean));  % 水平面目标，vz=0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. 添加环境白噪声（全频带）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_power = mean(p_clean.^2);  % 全频带信号功率
        noise_power = p_power / (10^(SNR/10));
        
        % 生成4通道独立白噪声（全频带）
        noise_p  = sqrt(noise_power) * randn(size(t));
        noise_vx = sqrt(noise_power) * randn(size(t));
        noise_vy = sqrt(noise_power) * randn(size(t));
        noise_vz = sqrt(noise_power) * randn(size(t));
        
        % 合成带噪信号（全频带）
        p_noisy  = p_clean  + noise_p;
        vx_noisy = vx_clean + noise_vx;
        vy_noisy = vy_clean + noise_vy;
        vz_noisy = vz_clean + noise_vz;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. 【关键】整体通过系统FIR（模拟接收链路带宽限制）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(Fir)
            p_n  = filter(Fir, 1, p_noisy);
            vx_n = filter(Fir, 1, vx_noisy);
            vy_n = filter(Fir, 1, vy_noisy);
            vz_n = filter(Fir, 1, vz_noisy);
        else
            p_n  = p_noisy;
            vx_n = vx_noisy;
            vy_n = vy_noisy;
            vz_n = vz_noisy;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5. 存储数据
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        idx_start = 1 + N*(iiii-1);
        idx_end   = iiii * N;
        Data(1, idx_start:idx_end) = vx_n;
        Data(2, idx_start:idx_end) = vy_n;
        Data(3, idx_start:idx_end) = vz_n;
        Data(4, idx_start:idx_end) = p_n;
    end
    
    Data = Data';  % 转置为 [时间 × 4] 格式
end