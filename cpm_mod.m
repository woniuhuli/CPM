function cpm_sig = cpm_mod(inputs,h,sps,L,q,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���CPM����
% input��������Ԫ
% h������ָ��
% sps��ÿ�����Ų�������
% L�����������䣩����
% q��������ֺ�������������ʽ������sps*L���㣩
% m�����ƽ�����
% cpm_sig��CPM�����ź�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_length = length(inputs); %���������Ŀ
if(L>1)
    inputs = [ones(1,L-1)*(1-m),inputs];  
    %���������0ʱ��֮ǰL-1������״̬��Ϊʲô��1-m?
end

theta_n = zeros(1,info_length+1);             %n=0�ĳ�ʼ��λ״̬
phases = zeros(1,sps*info_length);            %�ź���λֵ
for n = 0:info_length-1
    %�����n�����Ŷ�Ӧ����λֵҲ��[nT:(n+1)T]ʱ���
    In = n+L;%�������õ�n�����ŵ�Inputs�±�
    %t = [n*sps+1:(n+1)*sps];%��n������ռ�ݵ�ʱ��
    if n>0
        %��n�����Ŷ�Ӧ����λ״̬
        theta_n(n+1) = theta_n(n) + pi*h*inputs(In-L);  %In-L+1?
    end
    theta = zeros(1,sps);%����״̬
    for k = n-L+1:n
        %�����n�����ŵĹ���״̬
        theta = theta + 2*pi*h*inputs(k+L)*q((n-k)*sps+1:(n+1-k)*sps);
    end
    
    %������λ״̬�͹���״̬������λֵ
    phases(n*sps+1:(n+1)*sps) = mod(theta_n(n+1) + theta,2*pi);       
end
cpm_sig = exp(1j*phases);%������λֵ���㸴����
