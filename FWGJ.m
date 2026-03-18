clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设置参数 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_total = 50*60;           % 总时间 3分钟
fs = 512;                 % 采样频率
T = 10;                    % 每段时长 (秒)
N = fs * T;               % 每段点数
Num_N = T_total / T;      % 快拍数

t = (0:N-1) / fs;         % 时间向量 (长度=N)

frequency = [40, 50, 60, 70];  % 线谱频率
theta = [50, 200];             % 方位角 (度)
SNR = -10;                      % 信噪比 (dB)
deltaSigma = 1;                % 方位分辨率

p_sx = [0.45, 0.4, 0.4, 0.45]; % 线谱幅度

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 设计FIR滤波器（15–200 Hz带通）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fl = 15; fh = 200;
FIRLen = 256;
Fir = fir1(FIRLen, [fl/(fs/2), fh/(fs/2)], hamming(FIRLen+1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 生成仿真数据（✅ 关键：传入Fir！）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = Aim_Signal(N, Num_N, t, frequency, theta, SNR, p_sx, Fir);

% 提取四通道
ax_t = Data(:,1);   % vx
ay_t = Data(:,2);   % vy
az_t = Data(:,3);   % vz
p_t  = Data(:,4);   % 声压

% 截取第一段用于频谱观察
sim_x = Data(1:N, 1);
sim_y = Data(1:N, 2);
sim_p = Data(1:N, 4);

% FFT 频谱（单边谱）
Yp = fft(sim_p); P2p = abs(Yp / N); sim_p_fftp = 2*P2p(1:N/2);
Yx = fft(sim_x); P2x = abs(Yx / N); sim_x_fft  = 2*P2x(1:N/2);
Yy = fft(sim_y); P2y = abs(Yy / N); sim_y_fft  = 2*P2y(1:N/2);

f0 = (1:N/2) * (fs / N);  % 频率轴

% 绘图（可选）
figure('Name', '频域信号');
subplot(3,1,1); plot(f0, sim_p_fftp, 'k'); xlim([0 fs/2]); title('P');
subplot(3,1,2); plot(f0, sim_x_fft,  'k'); xlim([0 fs/2]); title('Vx');
subplot(3,1,3); plot(f0, sim_y_fft,  'k'); xlim([0 fs/2]); title('Vy');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 复声强方位估计
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_Data = p_t;
Vx_Data = ax_t;
Vy_Data = ay_t;

% Welch 参数
win_len = fs * 1;       % 1秒窗
overlap = 0.5;
step_len = floor(win_len * (1 - overlap));
Nfft = win_len;
f_fft = (0:Nfft-1) * fs / Nfft;
f_idx = find(f_fft >= fl & f_fft <= fh);
num_f = length(f_idx);

Old_PFFT  = zeros(Num_N, num_f);
Old_VXFFT = zeros(Num_N, num_f);
Old_VYFFT = zeros(Num_N, num_f);
Sigma_TF  = zeros(Num_N, num_f);

for iiii = 1:Num_N
    s1 = 1 + (iiii-1)*N;
    e1 = s1 + N - 1;
    
    ptemp  = P_Data(s1:e1);
    axtemp = Vx_Data(s1:e1);
    aytemp = Vy_Data(s1:e1);
    
    % === 直接使用已滤波信号（无需再滤波！）===
    Pf  = ptemp;
    Vxf = axtemp;
    Vyf = aytemp;
    
    % Welch 平均
    total_len = length(Pf);
    num_windows = floor((total_len - win_len) / step_len) + 1;
    
    Ptk = zeros(win_len, 1);
    Vxtk = zeros(win_len, 1);
    Vytk = zeros(win_len, 1);
    
    win = hanning(win_len);
    %%%welch操作。
    for xh = 1:num_windows
        start_idx = (xh-1)*step_len + 1;
        end_idx   = start_idx + win_len - 1;
        
        p_win  = Pf(start_idx:end_idx) .* win;
        vx_win = Vxf(start_idx:end_idx) .* win;
        vy_win = Vyf(start_idx:end_idx) .* win;
        
        Ptk  = Ptk  + fft(p_win, win_len);
        Vxtk = Vxtk + fft(vx_win, win_len);
        Vytk = Vytk + fft(vy_win, win_len);
    end
    
    Ptk  = Ptk  / num_windows;
    Vxtk = Vxtk / num_windows;
    Vytk = Vytk / num_windows;
    
    % 提取15-200Hz
    f_current = (0:win_len-1) * fs / win_len;
    f_idx_current = find(f_current >= fl & f_current <= fh);
    
    Old_PFFT(iiii, :)  = Ptk(f_idx_current).';
    Old_VXFFT(iiii, :) = Vxtk(f_idx_current).';
    Old_VYFFT(iiii, :) = Vytk(f_idx_current).';
    
    % 方位估计
    P  = Old_PFFT(iiii, :);
    Vx = Old_VXFFT(iiii, :);
    Vy = Old_VYFFT(iiii, :);
    
    Spx = conj(Vx) .* P;
    Spy = conj(Vy) .* P;
    Ix = real(Spx);
    Iy = real(Spy);
    
    sigma = atan2(Iy, Ix) * 180/pi;
    sigma(sigma < 0) = sigma(sigma < 0) + 360;
    Sigma_TF(iiii, :) = sigma;



%%%%%%%%%状态离散化，将方位角值转化为状态值
%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 新增：方位角按每10度映射为状态值（0-10°→1，11-20°→2，…，351-360°→36）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 初始化全局状态矩阵（和Sigma_TF维度完全一致：快拍数×频率数）
% 注意：只在第一次循环时初始化，避免重复清空！
if iiii == 1
    State_TF = zeros(size(Sigma_TF));  % 全局存储所有快拍的状态值
end

% 2. 取出当前快拍的所有方位角
current_sigma = Sigma_TF(iiii, :);  % 一维数组：当前快拍×所有频率的方位角

% 3. 核心：每10度映射为一个状态值（数学计算，无需手动写36个区间）
% 逻辑：(角度值-1) ÷ 10 取整 + 1 → 保证11-20°→2，351-360°→36
current_state = floor((current_sigma - 1) / 10) + 1;

% 4. 边界修正：处理360度（360°应归为状态36，而非37）
current_state(current_state > 36) = 36;
% 额外保护：防止角度为0的情况（0°归为状态1）
current_state(current_state < 1) = 1;

% 5. 将当前快拍的状态值存入全局矩阵
State_TF(iiii, :) = current_state;
end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 专利级状态转移熵计算（按【每个频点单独计算时间维度的转移熵】）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. 基础参数定义
num_freq = size(State_TF, 2);  % 频率点数量（num_f）
num_states = 36;               % 状态总数（1~36）
epsilon = 1e-10;               % 防止log2(0)的极小值
alpha = 0.4; beta = 0.6;       % 双熵融合权重

%% 2. 初始化存储矩阵（每个频点对应一组熵值）
% 列：频率点；行：熵类型（1=转移熵Htrans，2=分布熵Hdist，3=融合熵Hfusion）
entropy_results = zeros(3, num_freq);
% 存储每个频点的状态权重（可选，专利文档用）
state_weights_per_freq = zeros(num_states, num_freq);

%% 3. 逐频点计算转移熵（核心：同一频点，遍历所有时间/快拍）
for f = 1:num_freq
    % --------------------------
    % 步骤1：提取当前频点的时间序列（所有快拍的状态值）
    % --------------------------
    state_seq = State_TF(:, f);  % 一维数组：[快拍1, 快拍2, ..., 快拍Num_N]
    % 剔除无效值（若有）：只保留1~36的有效状态
    state_seq = state_seq(state_seq >= 1 & state_seq <= num_states);
    % 空序列保护（跳过无有效状态的频点）
    if length(state_seq) < 2  % 至少需要2个状态才能计算转移
        entropy_results(:, f) = [0, 0, 0];
        continue;
    end
    
    % --------------------------
    % 步骤2：统计当前频点的状态出现次数 & 权重
    % --------------------------
    C = zeros(1, num_states);  % Ci：当前频点状态i的出现次数
    for s = 1:num_states
        C(s) = sum(state_seq == s);
    end
    Ctotal = sum(C);           % 当前频点总状态数
    w = C / Ctotal;            % 当前频点各状态的权重（出现频率）
    state_weights_per_freq(:, f) = w;  % 存储权重（专利用）
    
    % --------------------------
    % 步骤3：构建当前频点的状态转移次数矩阵
    % --------------------------
    trans_count = zeros(num_states, num_states);
    % 遍历时间序列，统计i→j的转移次数（仅当前频点的时间维度）
    for k = 1:length(state_seq)-1
        i = state_seq(k);      % 当前快拍的状态
        j = state_seq(k+1);    % 下一个快拍的状态
        trans_count(i, j) = trans_count(i, j) + 1;
    end
    
    % --------------------------
    % 步骤4：计算当前频点各状态的条件转移熵 H(i)
    % --------------------------
    H = zeros(1, num_states);  % H(i)：状态i的条件转移熵
    for i = 1:num_states
        if C(i) == 0  % 未出现的状态，转移熵为0
            H(i) = 0;
            continue;
        end
        % 条件概率 P(j|i) = 转移次数(i→j) / 状态i总出现次数
        P_cond = trans_count(i, :) / C(i);
        P_cond(P_cond == 0) = epsilon;  % 避免log2(0)
        % 条件转移熵公式
        H(i) = -sum(P_cond .* log2(P_cond));
    end
    
    % --------------------------
    % 步骤5：计算当前频点的全局加权转移熵 Htrans
    % --------------------------
    Htrans = sum(w .* H);
    
    % --------------------------
    % 步骤6：计算当前频点的状态分布熵 Hdist（专利补充）
    % --------------------------
    P_dist = C / Ctotal;
    P_dist(P_dist == 0) = epsilon;
    Hdist = -sum(P_dist .* log2(P_dist));
    
    % --------------------------
    % 步骤7：双熵融合（专利级特征）
    % --------------------------
    Hfusion = alpha * Htrans + beta * Hdist;
    
    % --------------------------
    % 存储当前频点的熵值结果
    % --------------------------
    entropy_results(:, f) = [Htrans, Hdist, Hfusion];
end

%% 4. 结果输出（专利文档可直接引用）
fprintf('===== 逐频点状态转移熵计算结果 =====\n');
% 频率轴（对应entropy_results的列）
freq_axis = f_current(f_idx_current);  % 15~200Hz的频率值
for f = 1:num_freq
    fprintf('频率 %.1f Hz：转移熵=%.4f，分布熵=%.4f，融合熵=%.4f\n', ...
        freq_axis(f), entropy_results(1,f), entropy_results(2,f), entropy_results(3,f));
end

%% 5. 可视化（可选，专利附图用）
figure('Name','逐频点转移熵分布');
plot(freq_axis, entropy_results(1,:), 'k-o', 'LineWidth',1.2, 'MarkerSize',4);
xlabel('频率 / Hz'); ylabel('加权状态转移熵');
title('各频点时间维度的状态转移熵分布');
grid on; xlim([fl, fh]);

figure('Name','逐频点融合熵分布');
plot(freq_axis, entropy_results(3,:), 'r-s', 'LineWidth',1.2, 'MarkerSize',4);
xlabel('频率 / Hz'); ylabel('双熵融合值');
title('各频点时间维度的双熵融合特征');
grid on; xlim([fl, fh]);

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 熵值加权直方图法
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. 基础参数定义
deltaSigma = 1;                % 方位分辨率（1度）
Nsigma = round(360 / deltaSigma);  % 方位bin数量（360个）
az_axis = (0.5 : 1 : Nsigma - 0.5) * deltaSigma;  % 方位bin中心
num_target = 2;                % 预设目标数量（和theta=[50,200]一致）
true_az = [50, 200];           % 真实目标方位（用于精度分析）
Num_N = size(Sigma_TF, 1);     % 总时间片数（快拍数）
select_freq_num = 20;          % 新增：选定的频点数量（20个）

%% 2. 预计算频点熵值权重 + 筛选双熵最小的几个频点（核心修改）
epsilon = 1e-10;
alpha = 20;  % 放大系数，越大差异越明显（5=轻度，20=重度）

% 步骤1：按融合熵（双熵）升序排序，选出最小的20个频点
[fusion_entropy_sorted, fusion_entropy_idx] = sort(entropy_results(3, :));  % 融合熵升序
select_freq_idx = fusion_entropy_idx(1:min(select_freq_num, num_freq));    % 前20个低熵频点索引
fprintf('筛选出融合熵最小的%d个频点，对应频率：', length(select_freq_idx));
fprintf('%.1f ', freq_axis(select_freq_idx));
fprintf('\n');

% 步骤2：仅计算选中频点的熵值权重
entropy_weights_all = exp(-alpha * entropy_results(1, :));  % 所有频点的原始权重
entropy_weights = entropy_weights_all(select_freq_idx);     % 仅选中20个频点的权重
entropy_weights = entropy_weights / sum(entropy_weights);   % 归一化（仅针对20个频点）
num_freq_selected = length(entropy_weights);                % 实际选中的频点数（防止总频点不足20）

%% 3. 初始化结果存储（专利级结构化存储）
time_slice_target_az = zeros(Num_N, num_target);  
time_slice_error = zeros(Num_N, num_target);      
% time_slice_hist = zeros(Num_N, Nsigma);          

%% 4. 逐时间片计算加权直方图&目标方位（仅遍历选中的20个频点）
fprintf('===== 开始逐时间片计算目标方位（仅用%d个低熵频点）=====\n', num_freq_selected);
for t = 1:Num_N  % 遍历每个时间片（快拍1~Num_N）
    current_hist = zeros(1, Nsigma);  
    
    % --------------------------
    % 核心修改：仅遍历选中的20个频点
    % --------------------------
    for f_idx = 1:num_freq_selected
        f = select_freq_idx(f_idx);  % 原始频点索引
        freq_weight = entropy_weights(f_idx);  % 选中频点的权重
        az_val = Sigma_TF(t, f);               % 该时间片-该频点的方位角
        state_val = State_TF(t, f);            % 该时间片-该频点的状态值
        
        % 提取声强能量（无则用1）
        if exist('all_I', 'var') && size(all_I,1)>=t && size(all_I,2)>=f
            I_val = all_I(t, f) + epsilon;
        else
            I_val = 1;
        end
        
        % 无效值过滤
        if isnan(az_val) || state_val < 1 || state_val > 36 || I_val <= 0
            continue;
        end
        
        % 映射到方位bin + 加权累加
        bin = floor(az_val / deltaSigma) + 1;
        bin = min(max(bin, 1), Nsigma);
        current_hist(bin) = current_hist(bin) + I_val * freq_weight;
    end
    
% --------------------------
% 峰值检测 + 方位排序（原有逻辑保留）

% 直方图归一化（避免幅值影响峰值检测）
current_hist = current_hist / max(current_hist + epsilon);
% 存储当前时间片直方图（可选）
% time_slice_hist(t, :) = current_hist;

% 峰值检测（按能量降序取前num_target个峰值）
[peak_vals, peak_idx] = findpeaks(current_hist, ...
    'SortStr','descend', 'NPeaks', num_target, 'MinPeakHeight', 0.1);  % 过滤低能量峰值

% ===== 核心修复1：过滤掉peak_idx中的非正整数/NaN =====
peak_idx = peak_idx(~isnan(peak_idx) & peak_idx >= 1 & peak_idx <= length(az_axis));
peak_vals = peak_vals(~isnan(peak_idx) & peak_idx >= 1 & peak_idx <= length(az_axis));

% ===== 核心修复2：峰值补全（用前一帧结果，避免NaN）=====
if length(peak_idx) < num_target
    % 初始化补全的峰值索引（优先用前一帧结果，无则用0）
    fill_idx = [];
    if t > 1  % 非第一帧，用前一帧的目标方位对应的索引
        prev_az = time_slice_target_az(t-1, :);
        fill_idx = round(prev_az - 0.5);  % 方位角转索引（az_axis是0.5,1.5,...359.5）
        fill_idx = fill_idx(1:(num_target - length(peak_idx)));
    else  % 第一帧，无历史则用默认值（比如50°和200°对应的索引）
        fill_idx = [round(50-0.5), round(200-0.5)];
        fill_idx = fill_idx(1:(num_target - length(peak_idx)));
    end
    % 合并峰值索引（有效峰值 + 补全峰值）
    peak_idx = [peak_idx, fill_idx];
    peak_vals = [peak_vals, zeros(1, length(fill_idx))];  % 补全的峰值能量设为0
end

% 转换为方位角（此时peak_idx都是有效正整数）
current_target_az = az_axis(peak_idx);


% --------------------------
% 新增：按方位角大小排序（核心修正，仅3行）
% --------------------------
if length(current_target_az) >= 2  % 有2个目标时才排序
    [sorted_az, sorted_idx] = sort(current_target_az);  % 按方位角升序排列
    current_target_az = sorted_az;  % 小方位=目标1，大方位=目标2
    peak_vals = peak_vals(sorted_idx);  % 能量值同步排序（可选，不影响误差）
end

time_slice_target_az(t, :) = current_target_az;




    % 误差计算（原有逻辑保留）
    % --------------------------
    for i = 1:num_target
        if ~isnan(current_target_az(i))
            error1 = abs(current_target_az(i) - true_az(i));
            error2 = 360 - error1;
            time_slice_error(t, i) = min(error1, error2);
        else
            time_slice_error(t, i) = nan;
        end
    end
    
    % 进度输出
    if mod(t, 10) == 0
        fprintf('已完成 %d/%d 个时间片计算\n', t, Num_N);
    end
end

%% 5. 结果分析 & 可视化（原有逻辑完全保留）
% --------------------------
% 5.1 时间片-目标方位角曲线
% --------------------------
figure('Name','每个时间片的目标方位估计结果（20个低熵频点）');
subplot(2,1,1);
plot(1:Num_N, time_slice_target_az(:,1), 'r-o', 'LineWidth',1, 'MarkerSize',3);
hold on;
plot(1:Num_N, true_az(1)*ones(Num_N,1), 'r--', 'LineWidth',1.5);
xlabel('时间片编号'); ylabel('方位角 (deg)');
title(['第1个目标方位估计（真实值：', num2str(true_az(1)), '度）']);
grid on; ylim([0, 360]);
%ylim([true_az(1)-20, true_az(1)+20]);

subplot(2,1,2);
plot(1:Num_N, time_slice_target_az(:,2), 'b-s', 'LineWidth',1, 'MarkerSize',3);
hold on;
plot(1:Num_N, true_az(2)*ones(Num_N,1), 'b--', 'LineWidth',1.5);
xlabel('时间片编号'); ylabel('方位角 (deg)');
title(['第2个目标方位估计（真实值：', num2str(true_az(2)), '度）']);
grid on; ylim([0, 360]);
%ylim([true_az(2)-20, true_az(2)+20]);

% --------------------------
% 5.2 方位估计误差统计
% --------------------------
error1_valid = time_slice_error(:,1);
error1_valid = error1_valid(~isnan(error1_valid));
error2_valid = time_slice_error(:,2);
error2_valid = error2_valid(~isnan(error2_valid));

mean_error1 = mean(error1_valid);
std_error1 = std(error1_valid);
mean_error2 = mean(error2_valid);
std_error2 = std(error2_valid);

fprintf('\n===== 所有时间片方位估计误差分析=====\n');
fprintf('第1个目标（真实%.1f度）：平均误差=%.2f度，标准差=%.2f度\n', ...
    true_az(1), mean_error1, std_error1);
fprintf('第2个目标（真实%.1f度）：平均误差=%.2f度，标准差=%.2f度\n', ...
    true_az(2), mean_error2, std_error2);

% 误差分布直方图
figure('Name','方位估计误差分布（20个低熵频点）');
subplot(1,2,1);
histogram(error1_valid, 10);
xlabel('估计误差 (deg)'); ylabel('时间片数量');
title(['第1个目标误差分布（均值=%.2f）', mean_error1]);
grid on;

subplot(1,2,2);
histogram(error2_valid, 10);
xlabel('估计误差 (deg)'); ylabel('时间片数量');
title(['第2个目标误差分布（均值=%.2f）', mean_error2]);
grid on;

% --------------------------
% 5.3 估计稳定性分析（专利级）
% --------------------------
cv1 = std_error1 / mean_error1;
cv2 = std_error2 / mean_error2;

fprintf('\n===== 估计稳定性分析=====\n');
fprintf('第1个目标变异系数=%.3f（<1表示稳定）\n', cv1);
fprintf('第2个目标变异系数=%.3f（<1表示稳定）\n', cv2);

valid_rate1 = length(error1_valid) / Num_N * 100;
valid_rate2 = length(error2_valid) / Num_N * 100;
fprintf('第1个目标有效估计率=%.1f%%，第2个目标有效估计率=%.1f%%\n', ...
    valid_rate1, valid_rate2);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%绘图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 绘图
[T_idx, F_idx] = meshgrid(1:Num_N, f_current(f_idx_current));
figure;
scatter3(T_idx(:), F_idx(:), Sigma_TF(:), 15, 'k', 'filled');
% scatter3(T_idx(:), F_idx(:), State_TF(:), 15, 'k', 'filled');
xlabel('时间段'); ylabel('频率 / Hz'); zlabel('方位角 / deg');
title('时间-频率-方位角分布图');



figure;
plot(f_current(f_idx_current), Sigma_TF(end,:), 'k.', 'MarkerSize', 10);
% plot(f_current(f_idx_current), State_TF(end,:), 'k.', 'MarkerSize', 10);
xlabel('频率/Hz'); ylabel('方位角'); title('最终方位估计');

%%
%直方图法
% 参数
Nsigma = round(360 / deltaSigma);  % 确保整数
deltaSigma = 360 / Nsigma;         % 反推精确步长（避免浮点误差）

% 初始化直方图
X_pre = zeros(1, Nsigma);

% 遍历每个频率点的方位估计
for i = 1:length(sigma)
    ang = sigma(i);  % 当前方位角 ∈ [0, 360)
    
    % 【关键】映射到 bin 索引（1 到 Nsigma）
    bin = floor(ang / deltaSigma) + 1;
    
    % 安全保护：防止 ang == 360（理论上不会，但防万一）
    if bin > Nsigma
        bin = Nsigma;
    end
    
    % 累加声强能量
    X_pre(bin) = X_pre(bin) + sqrt(Ix(i)^2 + Iy(i)^2);
end

% 生成方位轴（每个 bin 的中心）
az_axis = (0.5 : 1 : Nsigma - 0.5) * deltaSigma;  % [0.5, 1.5, ..., 359.5]

% 绘图
figure('Name', '方位能量分布');
plot(az_axis, X_pre, 'k', 'LineWidth', 1.2);
xlabel('方位角 (deg)');
ylabel('归一化声强能量');
title('复声强方位估计结果');

xlim([0 360]);