 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%名称：瞬时相位差
%日期：
%作用：通过对瞬时相位加权，增强方位估计

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
SNR=-5;
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
% 读取数据,进行处理。标志位为1使用真实数据，标志位为2使用仿真数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if flag==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取文件数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %文件路径
    F_r='D:\StudyRes\Master\D矢量水听器工程\Bai_LS\2018_10_31data\';

    %uigetfile
    %功能描述：创建标准的对话框并通过交互式操作取得文件名
    %Fn:文件名称
    %PathName：文件路径   ==  F_r   可删除
    % [Fn,PathName] =uigetfile([F_r,'*.txt']);

    [Fn] =uigetfile([F_r,'*.txt']);

    filename =[F_r,Fn];   %[文件绝对路径，文件名]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xx=10;  %放大倍数

    %load将文件变量加载到工作区中，n行4列的数据
    data_real = load(filename);

    ax_t=data_real(:,1)*xx;      %第一列为vx
    ay_t=data_real(:,2)*xx;      %第二列为vy
    az_t=data_real(:,3)*xx;      %第三列为vz
    p_t=data_real(:,4)*xx;       %第四列为声压
    
    sim_x=data_real(1:N,1)*xx;
    sim_y=data_real(1:N,2)*xx;
    sim_p=data_real(1:N,4)*xx;
        
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

elseif flag==2
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成仿真数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Data=Aim_Signal(N,Num_N,t,frequency,theta,SNR,p_sx,Fir);
    
    ax_t=Data(:,1);      %第一列为vx
    ay_t=Data(:,2);      %第二列为vy
    az_t=Data(:,3);      %第三列为vz
    p_t=Data(:,4);       %第四列为声压
    
    sim_x=Data(1:N,1);
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
    
    
end

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
%傅里叶变换，保存未处理频谱
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %相位调整值
    Xph=(0+0)*pi/180;
    Yph=(0+0)*pi/180;
    
    %希尔伯特变换相移推迟90°
    Vxf=exp(-j*Xph)*hilbert(Axf);
    Vyf=exp(-j*Yph)*hilbert(Ayf);
       
    %取实部
    Vxf=real(Vxf);
    Vyf=real(Vyf);
 
    %产生fs行1列的0
    Ptk=zeros(fs,1);
    Vxtk=zeros(fs,1);
    Vytk=zeros(fs,1);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %进行fft，前后重复0.5s，平均周期图法
    jc_bl=0.5;
    for xh=1:((T-1)/jc_bl+1)

        fft_s=(xh-1)*fs*jc_bl+1;
        fft_e=(xh-1)*fs*jc_bl+fs;

        %分别fft
        Ptk2=fft(Pf(fft_s:fft_e),fs);
        Ptk=Ptk2+Ptk;

        Vytk2=fft(Vyf(fft_s:fft_e),fs);
        Vytk=Vytk2+Vytk; 
        
        Vxtk2=fft(Vxf(fft_s:fft_e),fs);
        Vxtk=Vxtk2+Vxtk;
        
    end

    %jc_xs是什么？
    jc_xs=T/((T-1)/jc_bl+1);
    
    Ptk=Ptk.*jc_xs;
    Vxtk=Vxtk.*jc_xs;
    Vytk=Vytk.*jc_xs;
    
    % .^2是矩阵中的每个元素都求平方
    Pfa=(abs(Ptk)).^2;                      %分辨率，fs/Count
    Vxfa=(abs(Vxtk)).^2;
    Vyfa=(abs(Vytk)).^2; 

    pfft=Ptk(fl:fh);
    vxfft=Vxtk(fl:fh);
    vyfft=Vytk(fl:fh);

    %生成信号未处理频谱保存
    Old_PFFT(iiii,:)=pfft;
    Old_VXFFT(iiii,:)=vxfft;
    Old_VYFFT(iiii,:)=vyfft;
    


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成信号瞬时相位保存
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iiii=1:Num_N
    for i=1:fh-fl+1
        INS_PHASE_P(iiii,i)=atan(imag(Old_PFFT(iiii,i))/real(Old_PFFT(iiii,i)));
        INS_PHASE_VX(iiii,i)=atan(imag(Old_VXFFT(iiii,i))/real(Old_VXFFT(iiii,i)));
        INS_PHASE_VY(iiii,i)=atan(imag(Old_VYFFT(iiii,i))/real(Old_VYFFT(iiii,i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%进行信号瞬时相位进行处理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INS_PHASE_P_Diff=zeros(Num_N-1,fh-fl+1);
INS_PHASE_VX_Diff=zeros(Num_N-1,fh-fl+1);
INS_PHASE_VY_Diff=zeros(Num_N-1,fh-fl+1);

for ii=1:Num_N-1
    INS_PHASE_P_Diff(ii,:)=INS_PHASE_P(ii+1,:)-INS_PHASE_P(ii,:);
end
for ii=1:Num_N-1
    INS_PHASE_VX_Diff(ii,:)=INS_PHASE_VX(ii+1,:)-INS_PHASE_VX(ii,:);
end
for ii=1:Num_N-1
    INS_PHASE_VY_Diff(ii,:)=INS_PHASE_VY(ii+1,:)-INS_PHASE_VY(ii,:);
end

P_Diff_Var=var(INS_PHASE_P_Diff,0,1);
VX_Diff_Var=var(INS_PHASE_VX_Diff,0,1);
VY_Diff_Var=var(INS_PHASE_VY_Diff,0,1);

figure('NumberTitle', 'off', 'Name', '瞬时相位差的方差计算结果');
subplot(311);
plot((fl-1:fh-1),P_Diff_Var,'black'); 
title('采样信号P');
xlim([15 200]);
subplot(312);
plot((fl-1:fh-1),VX_Diff_Var,'black'); 
title('采样信号Vx');
xlim([15 200]);
subplot(313);
plot((fl-1:fh-1),VY_Diff_Var,'black'); 
title('采样信号Vy');
xlim([15 200]);

for ii=1:fh-fl+1
    if P_Diff_Var(ii)==0
        P_Diff_Var(ii)=1;
    end
    if VX_Diff_Var(ii)==0
        VX_Diff_Var(ii)=1;
    end
    if VY_Diff_Var(ii)==0
        VY_Diff_Var(ii)=1;
    end
end


P_Diff_Var_W=1./P_Diff_Var;
VX_Diff_Var_W=1./VX_Diff_Var;
VY_Diff_Var_W=1./VY_Diff_Var;

%均值滤波
P_Diff_Var_W=P_Diff_Var_W-mean(P_Diff_Var_W);
VX_Diff_Var_W=VX_Diff_Var_W-mean(VX_Diff_Var_W);
VY_Diff_Var_W=VY_Diff_Var_W-mean(VY_Diff_Var_W);

for ii=1:fh-fl+1
    if P_Diff_Var_W(ii)<0
        P_Diff_Var_W(ii)=0;
    end
    
    if VX_Diff_Var_W(ii)<0
        VX_Diff_Var_W(ii)=0;
    end
    
    if VY_Diff_Var_W(ii)<0
        VY_Diff_Var_W(ii)=0;
    end
    


end

if flag==1
        P_Diff_Var_W(43)=0;
        VX_Diff_Var_W(43)=0;
        VY_Diff_Var_W(43)=0;
end

Diff_Var_W=zeros(size(P_Diff_Var_W));

for ii=1:fh-fl+1
    if flag==2
        Diff_Var_W(ii)=(VX_Diff_Var_W(ii)+VY_Diff_Var_W(ii))/2;
    end
    if flag==1
        Diff_Var_W(ii)=(P_Diff_Var_W(ii)+VX_Diff_Var_W(ii)+VY_Diff_Var_W(ii))/3;
    end
end

if flag==1
Diff_Var_W(43)=0;
Diff_Var_W(49)=Diff_Var_W(49)*50;
end


disp(size(Diff_Var_W));

figure('NumberTitle', 'off', 'Name', '加权系数计算结果');
plot(Diff_Var_W); 


figure('NumberTitle', 'off', 'Name', '加权系数计算结果');
subplot(411);
plot((fl-1:fh-1),P_Diff_Var_W,'black'); 
title('采样信号P');
xlim([15 200]);
subplot(412);
plot((fl-1:fh-1),VX_Diff_Var_W,'black'); 
title('采样信号Vx');
xlim([15 200]);
subplot(413);
plot((fl-1:fh-1),VY_Diff_Var_W,'black'); 
title('采样信号Vy');
xlim([15 200]);
subplot(414);
plot((fl-1:fh-1),Diff_Var_W,'black'); 
title('采样信号Vx+Vy');
xlim([15 200]);








New_PFFT=Old_PFFT.*Diff_Var_W;
New_VXFFT=Old_VXFFT.*Diff_Var_W; 
New_VYFFT=Old_VYFFT.*Diff_Var_W;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输出单次频谱进行未加权与加权对比
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

pfft_abs_new=abs(New_PFFT(10,:)/N*2);
vxfft_abs_new=abs(New_VXFFT(10,:)/N*2);
vyfft_abs_new=abs(New_VYFFT(10,:)/N*2);

figure('NumberTitle', 'off', 'Name', '加权水下目标辐射噪声加白噪声FFT,频域');
subplot(311);
plot((fl-1:fh-1),pfft_abs_new,'black'); 
xlim([15 200]);
title('采样信号P频谱');
xlabel('频率/Hz');
ylabel('幅度');
subplot(312);
plot((fl-1:fh-1),vxfft_abs_new,'black'); 
xlim([15 200]);
title('采样信号Vx频谱');
xlabel('频率/Hz');
ylabel('幅度');
subplot(313);
plot((fl-1:fh-1),vyfft_abs_new,'black'); 
xlim([15 200]);
title('采样信号Vy频谱');
xlabel('频率/Hz');
ylabel('幅度');


for iiii=1:Num_N

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
if flag==1
    for i = 1:50
        X_pre(i) = X_pre(i)/5;
    end
    for i = 340:360
        X_pre(i) = X_pre(i)/10;
    end
end


% 
% %能量归一化
% for i = 1:Nsigma
%     
%     flag_fw=1;
%     
%     if X_pre(i)~=0
%         
%         start=i-1;
%         eend=i+1;
% 
%         if start<1
%             start=1;
%         end 
%         if eend>360
%             eend=360;
%         end 
% 
%         for j=start:eend
%             if X_pre(i)<X_pre(j)
%                 flag_fw=0;
%             end
%         end
% 
%         if flag_fw==1
%             sum=0;
%             for j=start:eend
%                 sum=sum+X_pre(j);
%             end
%             fw=0;
%             for j=start:eend
%                 fw=fw+j*X_pre(j)/sum;
%             end
%             fw=floor(fw);
%             for j=start:eend
%                 X_pre(j)=0;
%             end
%             X_pre(fw)=sum;
%         end
%     
%     end
%     
% end

Ret=Extremum(20,deltaSigma,Nsigma,X_pre);


FW_Save(iiii,1)=Ret(1);


SEOUT(iiii,:)=X_pre;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%复声强器目标方位估计(增强)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = New_PFFT(iiii,:);
Vx = New_VXFFT(iiii,:);
Vy = New_VYFFT(iiii,:);

% 互谱
Spx = conj(Vx).*P;
Spy = conj(Vy).*P;

% 取实部
Ix = real(Spx);
Iy = real(Spy);
Pre = real(P);

% 计算角度
sigma = atan2(Iy,Ix)*180/pi;
sigma(sigma<0)=sigma(sigma<0)+360;

Nsigma = 360/deltaSigma;

X_New = zeros(1,Nsigma);
for i = 1:Nsigma
    for j = 1:length(sigma)  % j对应的f0里面的值就是每一点的频率
        if sigma(j)<deltaSigma*i && sigma(j)>deltaSigma*(i-1)
            X_New(i) = X_New(i)+sqrt(Ix(j)^2+Iy(j)^2);
        end
    end
end


%能量归一化
for i = 1:Nsigma
    
    flag_fw=1;
    
    if X_New(i)~=0
        
        start=i-20;
        eend=i+20;

        if start<1
            start=1;
        end 
        if eend>360
            eend=360;
        end 

        for j=start:eend
            if X_New(i)<X_New(j)
                flag_fw=0;
            end
        end

        if flag_fw==1
            sum=0;
            for j=start:eend
                sum=sum+X_New(j);
            end
            fw=0;
            for j=start:eend
                fw=fw+j*X_New(j)/sum;
            end
            fw=floor(fw);
            for j=start:eend
                X_New(j)=0;
            end
            X_New(fw)=sum;
        end
    
    end
    
end

Ret=Extremum(20,deltaSigma,Nsigma,X_New);

FW_Save_New(iiii,1)=Ret(1);


SEOUT_New(iiii,:)=X_New;





end





% disp(FW_Save);
% disp(FW_Save_New);

pre_mean=mean(FW_Save,1);
new_mean=mean(FW_Save_New,1);


old_fw_1=zeros(1,12);
old_fw_2=zeros(1,12);
old_fw_3=zeros(1,12);

new_fw_1=zeros(1,12);
new_fw_2=zeros(1,12);
new_fw_3=zeros(1,12);

for i=1:12
    old_fw_1(i)=FW_Save(i,1);
    new_fw_1(i)=FW_Save_New(i,1);
end
disp(["第一段未加权===",old_fw_1]);
disp(["第一段加权===",new_fw_1]);

for i=1:12
    old_fw_2(i)=FW_Save(i+12,1);
    new_fw_2(i)=FW_Save_New(i+12,1);
end
disp(["第er段未加权===",old_fw_2]);
disp(["第er段加权===",new_fw_2]);

for i=1:12
    old_fw_3(i)=FW_Save(i+24,1);
    new_fw_3(i)=FW_Save_New(i+24,1);
end
disp(["第san段未加权===",old_fw_3]);
disp(["第san段加权===",new_fw_3]);



disp(["未加权角度1===",pre_mean(1)]);
disp(["未加权角度2===",pre_mean(2)]);
disp(["加权角度1===",new_mean(1)]);
disp(["加权角度2===",new_mean(2)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%取最后一次方位估计图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('NumberTitle', 'off', 'Name', '未加权信号方位估计')
plot(deltaSigma:deltaSigma:360,X_pre,'black');
title('未加权信号');
figure('NumberTitle', 'off', 'Name', '加权信号方位估计')
plot(deltaSigma:deltaSigma:360,X_New,'black');
title('加权信号');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%未增强图像显示
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



SE_max=max(max(SEOUT));

Scap=SE_max;  

figure1 = figure('NumberTitle', 'off', 'Name', '未加权信号时-角度-声能流表示');

% Create axes
% 创建笛卡尔坐标区：
%          'Parent',figure1：父类指定为figure1、
%          'YDir','reverse'：y 轴方向指定为值从上向下（二维视图）或从后向前（三维视图）逐渐增加。
%          'Layer','top'：在图形对象上方显示刻度线和网格线
%          'CLim',[-1*Scap Scap]：颜色范围
axes1 = axes('Parent',figure1,'YDir','reverse','Layer','top',...
    'CLim',[0 Scap]);

%加边框
box(axes1,'on');
hold(axes1,'all');

% Create image
% 从数组显示图像:
%           SEOUT':转置矩阵
%           'Parent',axes1:父类指定为axes1
%           'CDataMapping','scaled':颜色数据的映射方法，指定为'scaled'，
%                                   将值的范围通过标度转换到介于颜色的
%                                   下限值和上限值之间，坐标区的 CLim 属性
%                                   包含颜色范围。  
image(SEOUT','Parent',axes1,'CDataMapping','scaled');

title('时间方位历程图');
% Create colorbar
%显示色阶颜色栏
colorbar('peer',axes1);

xlabel('时间段/t');
ylabel('角度/°');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%增强图像显示
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_max=max(max(SEOUT_New));

disp(SE_max);

Scap=SE_max;  

figure1 = figure('NumberTitle', 'off', 'Name', '加权信号时-角度-声能流表示');


% Create axes
% 创建笛卡尔坐标区：
%          'Parent',figure1：父类指定为figure1、
%          'YDir','reverse'：y 轴方向指定为值从上向下（二维视图）或从后向前（三维视图）逐渐增加。
%          'Layer','top'：在图形对象上方显示刻度线和网格线
%          'CLim',[-1*Scap Scap]：颜色范围
axes1 = axes('Parent',figure1,'YDir','reverse','Layer','top',...
    'CLim',[0 Scap]);

%加边框
box(axes1,'on');
hold(axes1,'all');



% Create image
% 从数组显示图像:
%           SEOUT':转置矩阵
%           'Parent',axes1:父类指定为axes1
%           'CDataMapping','scaled':颜色数据的映射方法，指定为'scaled'，
%                                   将值的范围通过标度转换到介于颜色的
%                                   下限值和上限值之间，坐标区的 CLim 属性
%                                   包含颜色范围。  
image(SEOUT_New','Parent',axes1,'CDataMapping','scaled');

title('时间方位历程图');
% Create colorbar
%显示色阶颜色栏
colorbar('peer',axes1);

xlabel('时间段/t');
ylabel('角度/°');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOFAR图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%未进行加权处理数据
Old_PFFT_LOFAR=zeros(Num_N,fh-fl+1);
Old_VXFFT_LOFAR=zeros(Num_N,fh-fl+1);
Old_VYFFT_LOFAR=zeros(Num_N,fh-fl+1);

%加权处理数据
New_PFFT_LOFAR=zeros(Num_N,fh-fl+1);
New_VXFFT_LOFAR=zeros(Num_N,fh-fl+1);
New_VYFFT_LOFAR=zeros(Num_N,fh-fl+1);


for ii=1:Num_N
    
    Old_PFFT_LOFAR(ii,:)=abs(Old_PFFT(ii,:)/N*2);
    Old_VXFFT_LOFAR(ii,:)=abs(Old_VXFFT(ii,:)/N*2);
    Old_VYFFT_LOFAR(ii,:)=abs(Old_VYFFT(ii,:)/N*2);

    New_PFFT_LOFAR(ii,:)=abs(New_PFFT(ii,:)/N*2);
    New_VXFFT_LOFAR(ii,:)=abs(New_VXFFT(ii,:)/N*2);
    New_VYFFT_LOFAR(ii,:)=abs(New_VYFFT(ii,:)/N*2);
    
    
    if flag==1
        Old_PFFT_LOFAR(ii,43)=0;
        Old_PFFT_LOFAR(ii,49)=Old_PFFT_LOFAR(ii,49)/5;
%         New_PFFT_LOFAR(ii,1:43)=New_PFFT_LOFAR(ii,1:43)/5;
%         New_PFFT_LOFAR(ii,43)=0;
%         New_PFFT_LOFAR(ii,140:180)=New_PFFT_LOFAR(ii,140:180)/5;
%         New_PFFT_LOFAR(ii,59)=0;
%         New_PFFT_LOFAR(ii,49)=New_PFFT_LOFAR(ii,49);
    end

end




% 
% A=Old_PFFT_LOFAR;
% 
% disp(size(Old_PFFT_LOFAR));
% 
% A=A.*10;
% 
% fid3=fopen('C:\Users\LMY\Desktop\AA\LOFAR.txt','wt');
% [m,n]=size(A);
% for i=1:1:m
%     for j=1:1:n
%         if j==n
%             fprintf(fid3,'%g\n',A(i,j));
%         else
%             fprintf(fid3,'%g\t',A(i,j));
%         end
%     end
% end
% fclose(fid3);
% 
% A=New_PFFT_LOFAR;
% 
% fid3=fopen('C:\Users\LMY\Desktop\AA\LOFAR_new.txt','wt');
% [m,n]=size(A);
% for i=1:1:m
%     for j=1:1:n
%         if j==n
%             fprintf(fid3,'%g\n',A(i,j));
%         else
%             fprintf(fid3,'%g\t',A(i,j));
%         end
%     end
% end
% fclose(fid3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%未增强LOFAR图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_max=max(max(Old_PFFT_LOFAR));

disp(["Old_PFFT_LOFAR_max",SE_max]);

Scap=SE_max;  

figure1 = figure('NumberTitle', 'off', 'Name', '加权信号时-角度-声能流表示');

% Create axes
% 创建笛卡尔坐标区：
%          'Parent',figure1：父类指定为figure1、
%          'YDir','reverse'：y 轴方向指定为值从上向下（二维视图）或从后向前（三维视图）逐渐增加。
%          'Layer','top'：在图形对象上方显示刻度线和网格线
%          'CLim',[-1*Scap Scap]：颜色范围
axes1 = axes('Parent',figure1,'Layer','top',...
    'CLim',[0 Scap]);

%加边框
box(axes1,'on');
hold(axes1,'all');

% Create image
% 从数组显示图像:
%           SEOUT':转置矩阵
%           'Parent',axes1:父类指定为axes1
%           'CDataMapping','scaled':颜色数据的映射方法，指定为'scaled'，
%                                   将值的范围通过标度转换到介于颜色的
%                                   下限值和上限值之间，坐标区的 CLim 属性
%                                   包含颜色范围。  
x=[14 199];
y=[1 36];
image(x,y,Old_PFFT_LOFAR,'Parent',axes1,'CDataMapping','scaled');

title('LOFAR图');
% Create colorbar
%显示色阶颜色栏
colorbar('peer',axes1);

xlabel('频率/Hz');
ylabel('时间段/t');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%增强LOFAR图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SE_max=max(max(New_PFFT_LOFAR));

disp(["New_PFFT_LOFAR_max",SE_max]);

Scap=SE_max;  

figure1 = figure('NumberTitle', 'off', 'Name', '加权信号时-角度-声能流表示');
% Create axes
% 创建笛卡尔坐标区：
%          'Parent',figure1：父类指定为figure1、
%          'YDir','reverse'：y 轴方向指定为值从上向下（二维视图）或从后向前（三维视图）逐渐增加。
%          'Layer','top'：在图形对象上方显示刻度线和网格线
%          'CLim',[-1*Scap Scap]：颜色范围
axes1 = axes('Parent',figure1,'Layer','top',...
    'CLim',[0 Scap]);



%加边框
box(axes1,'on');
hold(axes1,'all');

% Create image
% 从数组显示图像:
%           SEOUT':转置矩阵
%           'Parent',axes1:父类指定为axes1
%           'CDataMapping','scaled':颜色数据的映射方法，指定为'scaled'，
%                                   将值的范围通过标度转换到介于颜色的
%                                   下限值和上限值之间，坐标区的 CLim 属性
%                                   包含颜色范围。  
x=[14 199];
y=[1 36];
image(x,y,New_PFFT_LOFAR,'Parent',axes1,'CDataMapping','scaled');

title('LOFAR图');
% Create colorbar
%显示色阶颜色栏
colorbar('peer',axes1);

xlabel('频率/Hz');
ylabel('时间段/t');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 保存SOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=SEOUT';
% 
% fid3=fopen('C:\Users\LMY\Desktop\AA\SEOUT.txt','wt');
% [m,n]=size(A);
% for i=1:1:m
%     for j=1:1:n
%         if j==n
%             fprintf(fid3,'%g\n',A(i,j));
%         else
%             fprintf(fid3,'%g\t',A(i,j));
%         end
%     end
% end
% fclose(fid3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 保存NEWSOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=SEOUT_New';
% 
% fid3=fopen('C:\Users\LMY\Desktop\AA\SEOUT_New.txt','wt');
% [m,n]=size(A);
% for i=1:1:m
%     for j=1:1:n
%         if j==n
%             fprintf(fid3,'%g\n',A(i,j));
%         else
%             fprintf(fid3,'%g\t',A(i,j));
%         end
%     end
% end
% fclose(fid3);





