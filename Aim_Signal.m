%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 日期：2022/3/6（修正版：2025）
% 函数名称：生成矢量水听器接收目标信号
% 功能：双目标 + 多线谱 + 有色背景 + 白噪声干扰
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Data] = Aim_Signal(N, Num_N, t, frequency, theta, SNR, p_sx, Fir)

    Data = zeros(4, N * Num_N);
    deread = pi / 180;
    
    % 预计算滤波器（-6dB 噪声模型）
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
                  0.25*cos(2*pi*178*t)  + 0.21*cos(2*pi*100*t) + ...
                  0.28*cos(2*pi*85*t)  + 0.24*cos(2*pi*16*t) ;
                  
        background1 = filter(b_color, a_color, randn(size(t)));
        p1 = p1_line + background1;  % 背景作为目标自身辐射噪声
        
        vx1 = p1 * cos(theta(1) * deread);
        vy1 = p1 * sin(theta(1) * deread);
        
        % --- 目标2 ---
        p2_line = p_sx(3)*cos(2*pi*frequency(3)*t) + ...
                  p_sx(4)*cos(2*pi*frequency(4)*t) + ...
                  0.18*cos(2*pi*120*t)  + 0.26*cos(2*pi*30*t) + ...
                  0.29*cos(2*pi*148*t)  + 0.27*cos(2*pi*108*t) ;
                  
               
        background2 = filter(b_color, a_color, randn(size(t)));
        p2 = p2_line + background2;
        
        vx2 = p2 * cos(theta(2) * deread);
        vy2 = p2 * sin(theta(2) * deread);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2. 合成总信号（无外部噪声）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_total = p1 + p2;
        vx_total = vx1 + vx2;
        vy_total = vy1 + vy2;
        vz_total = zeros(size(p_total));  % 若无俯仰角，可设为0；否则需建模
        
        % 可选：通过FIR限制带宽（如只关心 20–200 Hz）
        if ~isempty(Fir)
            p_total  = filter(Fir, 1, p_total);
            vx_total = filter(Fir, 1, vx_total);
            vy_total = filter(Fir, 1, vy_total);
            vz_total = filter(Fir, 1, vz_total);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3. 添加环境白噪声（按标准 SNR 定义）
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_power = mean(p_total.^2);
        noise_power = p_power / (10^(SNR/10));  % 关键：按功率比例
        
        noise_p  = sqrt(noise_power) * randn(size(t));
        noise_vx = sqrt(noise_power) * randn(size(t));  % 假设各通道噪声功率相同
        noise_vy = sqrt(noise_power) * randn(size(t));
        noise_vz = sqrt(noise_power) * randn(size(t));
        
        p_n  = p_total  + noise_p;
        vx_n = vx_total + noise_vx;
        vy_n = vy_total + noise_vy;
        vz_n = vz_total + noise_vz;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4. 存储数据
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