function [complex_qam_data]=QAM16(bitdata)
%===============================================
%16QAM�ĵ��ƺ���������16QAM�ĸ����ź�
%���룺bitdata ��ʾ�������ݵ�������
%�����complex_qam_data ��ʾ��������������
%���ߣ��Ž�    2015��4��6��
%===============================================
X1=reshape(bitdata,4,length(bitdata)/4)';%���������ݰ������ң����ϵ��µ�˳��浽X1�����С�X1������4�У�ÿһ����1/4����������
d=1;%min distance of symble 
for i=1:length(bitdata)/4;
    for j=1:4
        X1(i,j)=X1(i,j)*(2^(4-j));%��ÿһ��4����������������һ��4λ�Ķ�������
    end
        source(i,1)=1+sum(X1(i,:));%%source(:,1)��ÿһ��ΪX1�����ж�Ӧ�е�4λ����������ֵ+1��ֵ��1-16��
end
mapping=[-3*d 3*d;
    	   -d  3*d;
            d  3*d;
            3*d  3*d;
            -3*d  d;
            -d  d;
            d  d;
            3*d  d;
            -3*d  -d; 
            -d  -d; 
            d  -d;
            3*d  -d;
            -3*d  -3*d;
            -d  -3*d;
            d  -3*d;
            3*d  -3*d];%mapping��һ��16*2�ľ���ÿһ�ж�Ӧ��һ��������ʵ�����鲿
 for i=1:length(bitdata)/4
     qam_data(i,:)=mapping(source(i),:);%qam_data��ÿһ��������Ԫ�أ��ֱ��ʾ��Ӧ���ݵ�QAMӳ������ʵ�����鲿
 end
 complex_qam_data=complex(qam_data(:,1),qam_data(:,2));%��qam_dataÿ�е�����Ԫ������������õ�����������complex_qam_data����Ϊ���

end

