function [ rcosw ] = rCosWindow( beta, Ts )
%===============================================
%���������Ҵ�����
%���룺beta ����ϵ��  �� Ts ����ѭ��ǰ׺��OFDM���ŵĳ��ȣ�Ϊ��ż��
%�����rconw ��ʾ��������������
%���ߣ��Ž�    2015��4��6��
%===============================================
t = 0 : (1 + beta) * Ts;
rcosw = zeros(1, (1 + beta) * Ts);
for i = 1 : beta * Ts;
    rcosw(i) = 0.5 + 0.5 * cos (pi + t(i) * pi / (beta * Ts));    
end
rcosw(beta * Ts + 1 : Ts) = 1;
for j = Ts + 1 : (1 + beta ) * Ts + 1;
    rcosw(j - i) = 0.5 + 0.5 * cos((t(j) - Ts) * pi / (beta * Ts));
end
rcosw = rcosw'; %�任Ϊ������
end

