function [next_state,pre_state] = state_grid(L,m,h_p,h_m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���״̬ת�������Ե�ǰ״̬����ǰ�������һ״̬��ʾ
% L�����䳤��
% m�����ƽ�����
% h_p,h_m������ָ������
% next_state����״Ϊ[state_num,m]����ʾ��ǰ״̬Ϊstate_index,��ǰ����Ϊinput_indexʱ����һ��ʱ�̵�state_index
% pre_state����״Ϊ[state_num,m], ��ʾ��ǰ״̬Ϊstate_index��ǰһʱ�̵�m�ֿ���״̬
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = h_m/h_p; %����ָ��
if mod(h_m,2)==0
    phase_state_num = h_p;                  %���ܵ���λ״̬��
else
    phase_state_num = 2*h_p;
end
phase_gap = 2*pi/phase_state_num;%������λ״̬��ļ��
phase_states = (0:phase_state_num-1)*phase_gap;%���ܵ���λ״̬

state_num = phase_state_num*(m^(L-1)); %���ܵ�״̬��
inputs = (0:m-1)*2-m+1;   %���ܵ��������
current_input_num = m;    %��ǰ������ſ�����
next_state = zeros(state_num,current_input_num);
pre_state = zeros(state_num,current_input_num);

if(L==1)
    %��ǰ״̬Ϊ[��n],state_index = phase_index
    %��һ״̬Ϊ[��n+1],��n+1 = ��n+��*h*In
    for state = 1:state_num
        phase_index = state;
        theta_n = phase_states(phase_index);
        for input = 1:current_input_num
            % ���㵱ǰ״̬��ǰ�������һ״̬��˳��Ҳ�ͼ��������һ״̬��ǰһ״̬
            In = inputs(input);
            theta_n_plaus1 = mod(theta_n + pi*h*In,2*pi);%��n+1
            next_state(state,input) = theta_n_plaus1/phase_gap+1;
            for j = 1:m
                if(pre_state(theta_n_plaus1/phase_gap+1,j) ==0)
                   pre_state(theta_n_plaus1/phase_gap+1,j) = state;
                   %continue
                   break
                end
            end
        end
    end
elseif(L==2)
    % ״̬Ϊ[��n,In_1],state_index = phase_index*m + In_1_index
    % ��һ״̬Ϊ[��n+1,In],��n+1 = ��n+��*h*In_1
    for state = 1:state_num
        phase_index = floor((state-1)/m)+1;%��n�±�
        theta_n = phase_states(phase_index);
        In_1_index = state-(phase_index-1)*m; %In-1�±�
        In_1 = inputs(In_1_index);
        for input = 1:current_input_num
            %In = inputs(input);
            theta_n_plaus1 = mod(theta_n + pi*h*In_1,2*pi);%��n+1
            theta_n_plaus1_index = theta_n_plaus1/phase_gap+1;%��n+1���±�
            next_state(state,input) = int8((theta_n_plaus1_index-1)*m+input);
            for j = 1:m
                if(pre_state(next_state(state,input),j)==0)
                    pre_state(next_state(state,input),j)=state;
                    break
                end
            end
        end
    end
elseif(L==3)
    % ״̬Ϊ[��n,In_1,In_2],state_index = phase_index*m^2 + In_1_index*m + In_2_index
    % ��һ״̬Ϊ[��n+1,In,In-1],��n+1 = ��n+��*h*In_2
    for state = 1:state_num
        phase_index = floor((state-1)/m^2)+1;
        theta_n = phase_states(phase_index);
        In_1_index = floor((state-(phase_index-1)*m^2-1)/m)+1;
        In_2_index = state-(phase_index-1)*m^2-(In_1_index-1)*m;
        In_2 = inputs(In_2_index);
        for input = 1:current_input_num
            theta_n_plaus1 = mod(theta_n + pi*h*In_2,2*pi);%��n+1
            theta_n_plaus1_index = theta_n_plaus1/phase_gap+1;%��n+1���±�
            next_state(state,input) = (theta_n_plaus1_index-1)*m^2+(input-1)*m+In_1_index;
            for j = 1:m
                if(pre_state(next_state(state,input),j)==0)
                    pre_state(next_state(state,input),j)=state;
                    break
                end
            end
        end
    end
                
    
elseif(L==4)
    % ״̬Ϊ[��n,In_1,In_2,In-3],state_index = phase_index*m^3 + In_1_index*m^2 + In_2_index*m + In_2_index
    % ��һ״̬Ϊ[��n+1,In,In-1,In_2],��n+1 = ��n+��*h*In_3
    for state = 1:state_num
        phase_index = floor((state-1)/m^3)+1;
        theta_n = phase_states(phase_index);
        In_1_index = floor((state-(phase_index-1)*m^3-1)/m^2)+1;
        In_2_index = floor((state-(phase_index-1)*m^3-(In_1_index-1)*m^2-1)/m)+1;
        In_3_index = state-(phase_index-1)*m^3-(In_1_index-1)*m^2-(In_2_index-1)*m;
        In_3 = inputs(In_3_index);
        for input = 1:current_input_num
            theta_n_plaus1 = mod(theta_n + pi*h*In_3,2*pi);%��n+1
            theta_n_plaus1_index = theta_n_plaus1/phase_gap+1;%��n+1���±�
            next_state(state,input) = (theta_n_plaus1_index-1)*m^3+(input-1)*m^2+(In_1_index-1)*m+(In_2_index);
            for j = 1:m
                if(pre_state(next_state(state,input),j)==0)
                    pre_state(next_state(state,input),j)=state;
                    break
                end
            end
        end
    end
end
end
