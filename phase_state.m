function correlator = phase_state(sps,m,h_m,h_p,q,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������㸽��������Ҫ�ĸ�����ֵ
% sps��sample per symbol
% m�����ƵĽ�����
% h_m,h_p������ָ������
% q�����庯������������ʽ����������Ϊsps*L
% L�����䣨����������
% correlator�����������Ļ���ؼ�������ı��ظ�����ֵ
% ��״Ϊ(state_num,m,:)��ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ�ĸ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = h_m/h_p; %����ָ��
if mod(h_m,2)==0
    phase_state_num = h_p;                  %���ܵ���λ״̬��
else
    phase_state_num = 2*h_p;
end
phase_gap = 2*pi/phase_state_num;%������λ״̬��ļ��
phase_states = (0:phase_state_num-1)*phase_gap;%���ܵ���λ״̬

% history_inputs = zeros(1,L-1); %L-1�п��ܵ���ʷ��������(In-1,In-2,In-L+1)
% for i = 1:m^(L-1)
%     %������ʷ���п��ܵ����̫�ָࣨ�����𣩣��ʺ���д��ͨʽ�����ǲ��þ���Lֵ�������
% end

state_num = phase_state_num*(m^(L-1)); %���ܵ�״̬��
inputs = (0:m-1)*2-m+1;   %���ܵ��������
current_input_num = m;    %��ǰ������ſ�����
correlator = zeros(state_num,current_input_num,sps);%�����ÿ������µ�һ������ʱ����ڵĸ�����ֵ

if(L == 1)
    % ��ǰ״̬Ϊ[theta_n]
    for state = 1:state_num
        % state_index = phase_index
        phase_index = state;
        theta_n = phase_states(phase_index);%��λ״̬��n
        for input = 1:current_input_num
            %����[state,input]��Ӧ�ĸ�����
            In = inputs(input);%��ǰ������ԪIn
            theta = 2*pi*h*In*q;%����״̬
            phase = mod(theta_n + theta,2*pi);%������λ״̬������״̬�����λֵ
            correlator(state,input,:) = exp(1j*phase);
        end
    end
elseif(L==2)
    % ��ǰ״̬Ϊ[theta_n,In_1]
    for state = 1:state_num
        % state_index = phase_index*m + In_1_index
        phase_index = floor((state-1)/m)+1;%״̬����mȡ������λ״̬�±�
        theta_n = phase_states(phase_index);
        In_1_index = mod(state-1,m) + 1;%״̬����mȡ���In-1�±�
        In_1 = inputs(In_1_index);%(n-1)Tʱ������
        for input = 1:current_input_num
            %����[state,input]��Ӧ�ĸ�����
            In = inputs(input);%��ǰ������ԪIn
            theta = 2*pi*h*In_1*q(sps+1:2*sps)+2*pi*h*In*q(1:sps);%����״̬
            phase = mod(theta_n + theta,2*pi);%������λ״̬������״̬�����λֵ
            correlator(state,input,:) = exp(1j*phase);
        end
    end    
elseif(L==3)
    % ��ǰ״̬Ϊ[theta_n,In_1,In_2]
    for state = 1:state_num
        % state_index = phase_index*m^2 + In_1_index*m + In_2_index
        phase_index = floor((state-1)/m^2)+1;
        theta_n = phase_states(phase_index);%��ǰ��λ״̬
        In_1_index = floor((state-(phase_index-1)*m^2-1)/m)+1;
        In_1 = inputs(In_1_index);%(n-1)Tʱ������
        In_2_index = state-(phase_index-1)*m^2-(In_1_index-1)*m;
        In_2 = inputs(In_2_index);%(n-2)Tʱ�̵�����
        for input = 1:current_input_num
            %����[state,input]��Ӧ�ĸ�����
            In = inputs(input);%nTʱ������
            theta = 2*pi*h*(In_2*q(2*sps+1:3*sps) + In_1*q(sps+1:2*sps) + In*q(1:sps));%����״̬
            phase = mod(theta_n + theta,2*pi);%������λ״̬������״̬�����λֵ
            correlator(state,input,:) = exp(1j*phase);
        end
    end
elseif(L==4)
    % ��ǰ״̬Ϊ[theta_n,In_1,In_2��In-3]
    for state = 1:state_num
        % state_index = phase_index*m^3 + In_1_index*m^2 + In_2_index*m + In_2_index
        phase_index = floor((state-1)/m^3)+1;
        theta_n = phase_states(phase_index);%��ǰ��λ״̬
        In_1_index = floor((state-(phase_index-1)*m^3-1)/m^2)+1;
        In_1 = inputs(In_1_index);%(n-1)Tʱ������
        In_2_index = floor((state-(phase_index-1)*m^3-(In_1_index-1)*m^2-1)/m)+1;
        In_2 = inputs(In_2_index);%(n-2)Tʱ�̵�����
        In_3_index = state-(phase_index-1)*m^3-(In_1_index-1)*m^2-(In_2_index-1)*m;
        In_3 = inputs(In_3_index);%(n-3)Tʱ�̵�����
        for input = 1:current_input_num
            %����[state,input]��Ӧ�ĸ�����
            In = inputs(input);%nTʱ������
            theta = 2*pi*h*(In_3*q(3*sps+1:4*sps)+In_2*q(2*sps+1:3*sps)+In_1*q(sps+1:2*sps)+In*q(1:sps));
            phase = mod(theta_n + theta,2*pi);%������λ״̬������״̬�����λֵ
            correlator(state,input,:) = exp(1j*phase);
        end
    end
end
        
