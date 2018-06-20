function viterbi_symbol = viterbi_demod(noisy_sig,correlator,next_states,L,h_m,h_p,m,sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CPM��ά�ر��������㷨
% noisy_sig�����ջ������ź�
% correlator���������֧������Ҫ��״̬���ظ����磬��״Ϊ(state_num,m,:)��ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ�ĸ�����
% pre_states��״̬ת��������״Ϊ[state_num,m],��ʾ��ǰ״̬Ϊstate_index,ǰһʱ��״̬m�ֿ��ܵ�state_index
% next_states��״̬ת��������״Ϊ[state_num,m]����ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ����һ��ʱ�̵�state_index
% L�����䳤��
% h_m��h_p������ϵ������
% m�����ƽ���
% sps��sample per symbol
% viterbi_symbol��ά�ر��������㷨����õ��ķ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(h_m,2)==0
    phase_state_num = h_p;                  %���ܵ���λ״̬��
else
    phase_state_num = 2*h_p;
end

state_num = phase_state_num*(m^(L-1)); %���ܵ�״̬��
input_num = m;
inputs = (0:m-1)*2-m+1;   %���ܵ��������

symbol_num = length(noisy_sig)/sps;

Infty = -1*1e10;
branch_metric = ones(state_num,m,3)*Infty;%ÿ��ʱ����ڵķ�֧����,����ǰ��ʼ��Ϊ������
%����[nT,(n+1)T]ʱ��Σ����ܵ�״̬ת�����Ϊ[state_num,m]�ֿ���
%(state_index,m,1)�����״̬ת�Ƶķ�֧����ֵ��(state_index,m,2)���״̬ת�ƶ�Ӧ��ǰһ״̬��
%(state_index,m,3)�����״̬ת�ƶ�Ӧ������

path_metric = zeros(symbol_num,state_num,3);%ȫʱ����ڵ�·������
%(symbol_index,state_index,1)����ǵ���symbol_index�����ŵĵ�state_index��״̬���ۼ�·������ֵ
%(symbol_index,state_index,2)����ǵ���symbol_index�����ŵĵ�state_index��״̬��ǰһʱ��״̬�����ڻ���
%(symbol_index,state_index,1)����ǵ���symbol_index�����ŵĵ�state_index��״̬��ǰһʱ������

%% ���ȳ�ʼ��t = 0ʱ��Ϊȫ0������t = [0:T]ʱ��εķ�֧������·��������Ҳ����һ����
y = noisy_sig(1:sps);%��һ�����ŵĽ����ź�
for state = 1:state_num
    for input = 1:input_num
       In = inputs(input);%��ǰ����
       R = reshape(correlator(state,input,:),1,sps); %״̬ת�ƶ�Ӧ�ı��ظ�����
       next_state = next_states(state,input);%��һ״̬
       metric_value = real(sum(y.*conj(R)));%�뱾�ظ���������ػ��״̬ת�ƵĶ���ֵ
       for j = 1:m
           if(branch_metric(next_state,j,1) == Infty)
               branch_metric(next_state,j,1)= metric_value;%״̬ת�Ʒ�֧����ֵ
               branch_metric(next_state,j,2) = state;%ת̬ת�Ƶ�ǰһ״̬
               branch_metric(next_state,j,3) = In;%״̬ת�ƶ�Ӧ������
               break
           end
       end
    end
end
temp = reshape(branch_metric(:,:,1),state_num,m);%��ȡ��֧����ֵ
[max_metric,index] = max(temp,[],2);%�ҵ�������ֵ����Ӧ���±�
path_metric(1,:,1) = max_metric;%��һ��symbol��״̬�Ķ���ֵ
for i = 1:state_num
    path_metric(1,i,2) = branch_metric(i,index(i),2);%ת�Ƶ���һ��symbol�Ҵ��ǰһ��״̬
    path_metric(1,i,3) = branch_metric(i,index(i),3);%ת�Ƶ���һ��symbol�Ҵ������
end

%% ��ǰ�������ÿ��ʱ��εķ�֧����ֵ�Լ��ۼ�·������ֵ
for t = 2:symbol_num
    %�����t�����Ŷ�Ӧ�ķ�֧������·������
    branch_metric = ones(state_num,m,3)*Infty;  %��֧��������ǰ��ʼ��Ϊ������
    y = noisy_sig((t-1)*sps+1:t*sps);           %��t�����Ŷ�Ӧ������
    for state = 1:state_num
        for input = 1:input_num
            In = inputs(input);%��ǰ����
            R = reshape(correlator(state,input,:),1,sps); %״̬ת�ƶ�Ӧ�ı��ظ�����
            next_state = next_states(state,input);%��һ״̬
            metric_value = real(sum(y.*conj(R)));%�뱾�ظ���������ػ��״̬ת�Ƶķ�֧����ֵ
            for j = 1:m
               if(branch_metric(next_state,j,1) == Infty)
                   branch_metric(next_state,j,1)= metric_value + path_metric(t-1,state,1);%��ǰת�Ʒ�֧����ֵ+ǰһ״̬�ۼ�·������ֵ
                   branch_metric(next_state,j,2) = state;%ת̬ת�Ƶ�ǰһ״̬
                   branch_metric(next_state,j,3) = In;%״̬ת�ƶ�Ӧ������
                   break
               end
            end
        end
    end
    temp = reshape(branch_metric(:,:,1),state_num,m);%��ȡ��֧����ֵ
    [max_metric,index] = max(temp,[],2);%�ҵ�������ֵ����Ӧ���±�
    path_metric(t,:,1) = max_metric;%��t��symbol��״̬�Ķ���ֵ
    for i = 1:state_num
        path_metric(t,i,2) = branch_metric(i,index(i),2);%ת�Ƶ���t��symbol�Ҵ��ǰһ��״̬
        path_metric(t,i,3) = branch_metric(i,index(i),3);%ת�Ƶ���t��symbol�Ҵ������
    end
end

%% Ȼ�����·������ֵ�Ӻ���ǰ����
viterbi_symbol = zeros(1,symbol_num);
% �����ҵ����·�����յ�
temp = reshape(path_metric(symbol_num,:,1),1,state_num);%����ʱ�̸���״̬���ۼ�·������ֵ
[~,index] = max(temp);%��������ʱ�̸�״̬���ۼ�·������ֵ�õ�����ʱ�̵����״̬
viterbi_symbol(symbol_num) = path_metric(symbol_num,index,3);%����ʱ�����״̬��Ӧ��������Ϊ����ʱ������Ĺ���ֵ
pre_state = path_metric(symbol_num,index,2);%����ʱ�����״̬��ǰһ״̬���ڻ���
%Ȼ�����·��ͼ�Ӻ���ǰ����
for t = symbol_num-1:-1:1
    viterbi_symbol(t) = path_metric(t,pre_state,3); %����״̬�õ�ǰ�������
    pre_state = path_metric(t,pre_state,2);         %״̬����         
end













