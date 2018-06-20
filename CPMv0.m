%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPM调制解调流程
% 这里不考虑载波部分，仅仅研究其等效基带信号，也就是复包络（若不研究载波，两者可认为等价）
% 2018/6/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 参数设置
m = 4;          %进制数
L = 2;          %关联长度，记忆长度
h_m = 2;h_p = 3;
h = h_m/h_p;    %调制指数
sps = 8;        %每个符号样点数 sample per symbol
N = 256;        %原始码元个数

%% CPM调制
[g,q] = rc_pulse(sps,L);%生成升余弦脉冲函数g，及其积分函数q
info_bit = randi([0,1],1,N);   %信息码元
temp = reshape(info_bit,log2(m),N/log2(m));%串并转换
temp = temp';
temp = bi2de(temp,'left-msb');
symbol = (2*temp-m+1)';%码元符号,行向量
CPM_sig = cpm_mod(symbol,h,sps,L,q,m);

%% 过AWGN信道
EbN0 = 10;       %dB
SNR = EbN0  - 10 * log10(sps) + 10 * log10(log2(m));  %EbN0--->SNR  计算方式?如下(2018/6/19)
%SNR = 10log10(S/N),S = Eb*log2(m)/sps每个样点能量值，也即信号功率，N = N0噪声功率
noisy_sig = awgn(CPM_sig,SNR,'measured');
%noisy_sig = CPM_sig;

%% 维特比解调算法

% 首先得到计算度量值需要的pM^L(或者2pM^L)个可能的相位值的复包络
% I = (In,In-1,In-2...In-L+1)长为L个序列，每个值有m种可能，θ(n)有p（2p)种可能
correlator = phase_state(sps,m,h_m,h_p,q,L);
%形状为(state_num,m,:)表示当前状态为state_index,当前输入为input_index时的复包络

% 然后得到状态网络图，每个时刻的状态为[θn,In-1,In-2,...In-L+1]
% 当前可能的输入有m种情况，因此总的转移状态情况为[state_num,m]
[next_states,pre_states] = state_grid(L,m,h_p,h_m); 
% next_state形状为[state_num,m]，表示当前状态为state_index,当前输入为input_index时的下一个时刻的state_index
% pre_state形状为[state_num,m],表示当前状态为state_index,前一时刻状态m种可能的state_index

% 维特比译码算法
viterbi_symbol = viterbi_demod(noisy_sig,correlator,next_states,L,h_m,h_p,m,sps);%维特比译码符号
symbol_error_num =  length(find(viterbi_symbol ~= symbol));%误符号数
symbol_error_radio = symbol_error_num/length(symbol);%误符号率
temp = (viterbi_symbol-1+m)/2;
temp = de2bi(temp,'left-msb');
temp = temp';
viterbi_bit = reshape(temp,1,[]);%维特比译码解调所得bit
[bit_error_num,bit_error_radio] = biterr(info_bit,viterbi_bit);%误码数及误码率


