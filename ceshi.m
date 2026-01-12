close all;
warning off
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置各种参数 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%标志位flag（=1，使用真实数据，=2使用仿真数据）
flag=2;


T_total=3*60;           %采样总时间 3分钟*60s
fs=512;                 %采样频率   512  1024    2的倍数

T=5;                    %每段采样时间/s
N=fs*T;                 %每段采样点数

Num_N=T_total/T;        %总的采样段数目

f0 = fs/N*(0:N-1);      %fft变换后每一段频率

t=0:1/fs:T-1/fs;

%方便角度转换
deread=pi/180;
redge=180/pi;

%两个信号频率
frequency=[40,50,60,70];

%两个角的入射角度
theta=[50,200];
%信噪比。
SNR=-15;
%方位估计精度
deltaSigma=1;

% 第一个目标添加线谱系数
p1_sx1=0.5;
p1_sx2=0.5;
% 第二个目标添加线谱系数
p2_sx1=0.5;
p2_sx2=0.5;

p_sx=[p1_sx1,p1_sx2,p2_sx1,p2_sx2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成滤波器(生成目标信号滤波器)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRLen：滤波器的阶数
%[wn1，wn2]：截止频率为采样率的一半，例如：15~200hz，对应fl/(fs/2),fh/(fs/2)
%hamming(FIRLen+1)：窗函数的长度应等于FIR 滤波器系数个数，即阶数 n+1。
fl=15;fh=200;
FIRLen=256;
Fir=fir1(FIRLen,[fl/(fs/2),fh/(fs/2)],hamming(FIRLen+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成仿真数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Data=Aim_Signal(N,Num_N,t,frequency,theta,SNR,p_sx,Fir);
    
    ax_t=Data(:,1);      %第一列为vx
    ay_t=Data(:,2);      %第二列为vy
    az_t=Data(:,3);      %第三列为vz
    p_t=Data(:,4);       %第四列为声压
    
    sim_x=Data(1:N,1);%第一段傅里叶变换
    sim_y=Data(1:N,2);
    sim_p=Data(1:N,4);
        
    sim_p_fft=abs(fft(sim_p)/N*2);
    sim_x_fft=abs(fft(sim_x)/N*2);
    sim_y_fft=abs(fft(sim_y)/N*2);

    figure('NumberTitle', 'off', 'Name', '仿真信号生成-频域表示');
    subplot(311);
    plot(f0,sim_p_fft,'black'); 
    title('采样信号P');
    subplot(312);
    plot(f0,sim_x_fft,'black'); 
    title('采样信号Vx');
    subplot(313);
    plot(f0,sim_y_fft,'black'); 
    title('采样信号Vy');
    









    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成滤波器、姿态偏转角度等
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%矢量悬挂偏角
ZT_FW=0;
%%%%%%%%%%%%300~1000  
fl=15;fh=200;
FIRLen=256;    %fir滤波器阶数
%floor：取比它小的整数，即朝负无穷方向取整
FIRLenHalf=floor(FIRLen/2);
%生成滤波器
%FIRLen：滤波器的阶数
%[wn1，wn2]：截止频率为采样率的一半，例如：15~200hz，对应fl/(fs/2),fh/(fs/2)
%hamming(FIRLen+1)：窗函数的长度应等于FIR 滤波器系数个数，即阶数 n+1。
coeff=fir1(FIRLen,[fl/(fs/2),fh/(fs/2)],hamming(FIRLen+1));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%设置存储数据的数组
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%保存相位
INS_PHASE_P=zeros(Num_N,fh-fl+1);
INS_PHASE_VX=zeros(Num_N,fh-fl+1);
INS_PHASE_VY=zeros(Num_N,fh-fl+1);

%未进行加权处理数据
Old_PFFT=zeros(Num_N,fh-fl+1);
Old_VXFFT=zeros(Num_N,fh-fl+1);
Old_VYFFT=zeros(Num_N,fh-fl+1);

%加权处理数据
New_PFFT=zeros(Num_N,fh-fl+1);
New_VXFFT=zeros(Num_N,fh-fl+1);
New_VYFFT=zeros(Num_N,fh-fl+1);

%输出 “时-方位-声能流” 图
SEOUT=zeros(Num_N,360/deltaSigma);
SEOUT_New=zeros(Num_N,360/deltaSigma);

%存储角度
FW_Save=zeros(Num_N,2);
FW_Save_New=zeros(Num_N,2);







    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 开始处理数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%  
P_Data=p_t;
Vz_Data=az_t;    
%%%%%%%%%%%%%%%为什么取反？
Vx_Data=ax_t;   
Vy_Data=ay_t;    

for iiii=1:Num_N
    
    %计算出在整个data中每一段数据的开头和结尾 
    s1=1+(iiii-1)*N;
    e1=s1+N-1;
    
    %将每一段数据临时保存下来
    ptemp=P_Data(s1:e1);   
    axtemp=Vx_Data(s1:e1);
    aytemp=Vy_Data(s1:e1);
    
    %mean求出该列的平均值
    %均值滤波
%     ptemp=ptemp-mean(ptemp);
%     axtemp=axtemp-mean(axtemp);
%     aytemp=aytemp-mean(aytemp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 目标加滤波器（无噪声）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %conv与filter区别导致需要截取数据
    
    pf=conv(coeff,ptemp);
    axf=conv(coeff,axtemp);
    ayf=conv(coeff,aytemp);
    
    %分段数据的长度
    TempLen=length(pf);
    %disp(['分段数据长度=    ',num2str(TempLen)]);
    
    %截取滤波后的数据
    Pf=pf(FIRLenHalf+1:TempLen-FIRLenHalf);
    Axf=axf(FIRLenHalf+1:TempLen-FIRLenHalf);
    Ayf=ayf(FIRLenHalf+1:TempLen-FIRLenHalf);
    
    PfLen=length(Pf);
    %disp(['分段截取后数据长度=    ',num2str(PfLen)]);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 傅里叶变换（重叠平均周期图法，修正版）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 去掉多余的希尔伯特变换（直接用滤波后的信号）
Vxf = Axf;  % 滤波后的Vx信号（时域，实数）
Vyf = Ayf;  % 滤波后的Vy信号（时域，实数）
Pf = Pf;    % 滤波后的P信号（时域，实数）

% 2. 重叠平均参数设置（明确变量含义，避免索引错误）
fs = 512;          % 采样频率
win_len = fs * 1;  % 每段FFT的窗口长度（1秒，可调整，如0.5s=256点）
overlap = 0.5;     % 重叠比例（0~1，如0.5=50%重叠）
overlap_len = floor(win_len * overlap);  % 重叠长度（如256点）
step_len = win_len - overlap_len;        % 每段步长（如256点）

% 3. 计算循环次数（确保不超出信号长度）
total_len = length(Pf);  % 信号总长度=N=2560
num_windows = floor((total_len - win_len) / step_len) + 1;  % 有效窗口数（如(2560-512)/256 +1=9）

% 4. 初始化存储变量（复数数组，保存FFT结果）
Ptk = zeros(fs, 1);  % FFT点数=fs（或win_len，推荐win_len）
Vxtk = zeros(fs, 1);
Vytk = zeros(fs, 1);

% 5. 重叠平均FFT（循环不越界）
for xh = 1:num_windows
    % 计算当前窗口的索引（确保不超出信号长度）
    fft_s = (xh - 1) * step_len + 1;
    fft_e = fft_s + win_len - 1;  % 窗口长度=win_len（如512点）
    
    % 截取当前窗口数据（长度=win_len）
    p_win = Pf(fft_s:fft_e);
    vx_win = Vxf(fft_s:fft_e);
    vy_win = Vyf(fft_s:fft_e);
    
    % 加窗（可选，减少频谱泄漏，推荐hanning窗）
    win = hanning(win_len);  % 汉宁窗（长度=win_len）
    p_win = p_win .* win;
    vx_win = vx_win .* win;
    vy_win = vy_win .* win;
    
    % FFT（点数=win_len，与窗口长度一致，避免补0）
    p_fft = fft(p_win, win_len);
    vx_fft = fft(vx_win, win_len);
    vy_fft = fft(vy_win, win_len);
    
    % 叠加FFT结果
    Ptk(1:win_len) = Ptk(1:win_len) + p_fft;
    Vxtk(1:win_len) = Vxtk(1:win_len) + vx_fft;
    Vytk(1:win_len) = Vytk(1:win_len) + vy_fft;
end

% 6. 归一化（修正重叠平均的幅值）
Ptk = Ptk / num_windows;  % 除以窗口数，抵消叠加放大
Vxtk = Vxtk / num_windows;
Vytk = Vytk / num_windows;

% 7. 计算功率谱（可选，用于观察）
Pfa = abs(Ptk).^2;
Vxfa = abs(Vxtk).^2;
Vyfa = abs(Vytk).^2;

% 8. 提取15~200Hz的频点（用频率索引，避免硬编码）
f0 = fs / win_len * (0:win_len-1);  % 频率轴（与FFT点数匹配）
f_idx = find(f0 >= fl & f0 <= fh);  % 15~200Hz对应的索引
pfft = Ptk(f_idx);
vxfft = Vxtk(f_idx);
vyfft = Vytk(f_idx);

% 9. 保存频谱
Old_PFFT(iiii,:) = pfft;
Old_VXFFT(iiii,:) = vxfft;
Old_VYFFT(iiii,:) = vyfft;

end



pfft_abs=abs(Old_PFFT(10,:)/N*2);
vxfft_abs=abs(Old_VXFFT(10,:)/N*2);
vyfft_abs=abs(Old_VYFFT(10,:)/N*2);

figure('NumberTitle', 'off', 'Name', '未加权水下目标辐射噪声加白噪声FFT,频域');
subplot(311);
plot((fl-1:fh-1),pfft_abs,'black'); 
xlim([15 200]);
title('采样信号P频谱');
xlabel('频率/Hz');
ylabel('幅度');
subplot(312);
plot((fl-1:fh-1),vxfft_abs,'black'); 
xlim([15 200]);
title('采样信号Vx频谱');
xlabel('频率/Hz');
ylabel('幅度');
subplot(313);
plot((fl-1:fh-1),vyfft_abs,'black'); 
xlim([15 200]);
title('采样信号Vy频谱');
xlabel('频率/Hz');
ylabel('幅度');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%复声强器目标方位估计(未增强)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = Old_PFFT(iiii,:);
Vx = Old_VXFFT(iiii,:);
Vy = Old_VYFFT(iiii,:);
% 互谱
Spx = conj(Vx).*P;
Spy = conj(Vy).*P;

% 取实部
Ix = real(Spx);
Iy = real(Spy);
% Pre = real(P);

% 计算角度
sigma = atan2(Iy,Ix)*180/pi;
sigma(sigma<0)=sigma(sigma<0)+360;

Nsigma = 360/deltaSigma;

X_pre = zeros(1,Nsigma);
for i = 1:Nsigma
    for j = 1:length(sigma)  % j对应的f0里面的值就是每一点的频率
        if sigma(j)<deltaSigma*i && sigma(j)>deltaSigma*(i-1)
            X_pre(i) = X_pre(i)+sqrt(Ix(j)^2+Iy(j)^2);
        end
    end
end





figure('NumberTitle', 'off', 'Name', '未加权水下目标辐射噪声加白噪声FFT,频域');

plot(15 :1:200, sigma, 'k.', 'MarkerSize', 10); 

title('二维频率方位角图');
xlabel('频率/Hz');
ylabel('方位角');







for i = 1:Nsigma

    flag_fw=1;

    if X_pre(i)~=0

        start=i-1;
        eend=i+1;

        if start<1
            start=1;
        end 
        if eend>360
            eend=360;
        end 

        for j=start:eend
            if X_pre(i)<X_pre(j)
                flag_fw=0;
            end
        end

        if flag_fw==1
            sum=0;
            for j=start:eend
                sum=sum+X_pre(j);
            end
            fw=0;
            for j=start:eend
                fw=fw+j*X_pre(j)/sum;
            end
            fw=floor(fw);
            for j=start:eend
                X_pre(j)=0;
            end
            X_pre(fw)=sum;
        end

    end

end



figure('NumberTitle', 'off', 'Name', '删除频点信号方位估计')
plot(deltaSigma:deltaSigma:360,X_pre,'black');
title('删除频点信号');