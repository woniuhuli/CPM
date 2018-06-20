%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPM���ƽ������
% ���ﲻ�����ز����֣������о����Ч�����źţ�Ҳ���Ǹ����磨�����о��ز������߿���Ϊ�ȼۣ�
% 2018/6/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ��������
m = 4;          %������
L = 2;          %�������ȣ����䳤��
h_m = 2;h_p = 3;
h = h_m/h_p;    %����ָ��
sps = 8;        %ÿ������������ sample per symbol
N = 256;        %ԭʼ��Ԫ����

%% CPM����
[g,q] = rc_pulse(sps,L);%�������������庯��g��������ֺ���q
info_bit = randi([0,1],1,N);   %��Ϣ��Ԫ
temp = reshape(info_bit,log2(m),N/log2(m));%����ת��
temp = temp';
temp = bi2de(temp,'left-msb');
symbol = (2*temp-m+1)';%��Ԫ����,������
CPM_sig = cpm_mod(symbol,h,sps,L,q,m);

%% ��AWGN�ŵ�
EbN0 = 10;       %dB
SNR = EbN0  - 10 * log10(sps) + 10 * log10(log2(m));  %EbN0--->SNR  ���㷽ʽ?����(2018/6/19)
%SNR = 10log10(S/N),S = Eb*log2(m)/spsÿ����������ֵ��Ҳ���źŹ��ʣ�N = N0��������
noisy_sig = awgn(CPM_sig,SNR,'measured');
%noisy_sig = CPM_sig;

%% ά�رȽ���㷨

% ���ȵõ��������ֵ��Ҫ��pM^L(����2pM^L)�����ܵ���λֵ�ĸ�����
% I = (In,In-1,In-2...In-L+1)��ΪL�����У�ÿ��ֵ��m�ֿ��ܣ���(n)��p��2p)�ֿ���
correlator = phase_state(sps,m,h_m,h_p,q,L);
%��״Ϊ(state_num,m,:)��ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ�ĸ�����

% Ȼ��õ�״̬����ͼ��ÿ��ʱ�̵�״̬Ϊ[��n,In-1,In-2,...In-L+1]
% ��ǰ���ܵ�������m�����������ܵ�ת��״̬���Ϊ[state_num,m]
[next_states,pre_states] = state_grid(L,m,h_p,h_m); 
% next_state��״Ϊ[state_num,m]����ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ����һ��ʱ�̵�state_index
% pre_state��״Ϊ[state_num,m],��ʾ��ǰ״̬Ϊstate_index,ǰһʱ��״̬m�ֿ��ܵ�state_index

% ά�ر������㷨
viterbi_symbol = viterbi_demod(noisy_sig,correlator,next_states,L,h_m,h_p,m,sps);%ά�ر��������
symbol_error_num =  length(find(viterbi_symbol ~= symbol));%�������
symbol_error_radio = symbol_error_num/length(symbol);%�������
temp = (viterbi_symbol-1+m)/2;
temp = de2bi(temp,'left-msb');
temp = temp';
viterbi_bit = reshape(temp,1,[]);%ά�ر�����������bit
[bit_error_num,bit_error_radio] = biterr(info_bit,viterbi_bit);%��������������


