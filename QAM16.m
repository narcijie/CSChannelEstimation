function [complex_qam_data]=QAM16(bitdata)
%===============================================
%16QAM的调制函数。产生16QAM的复数信号
%输入：bitdata 表示输入数据的行向量
%输出：complex_qam_data 表示输出结果的列向量
%作者：杜捷    2015年4月6日
%===============================================
X1=reshape(bitdata,4,length(bitdata)/4)';%将输入数据按从左到右，从上到下的顺序存到X1矩阵中。X1矩阵有4列，每一列有1/4个输入数据
d=1;%min distance of symble 
for i=1:length(bitdata)/4;
    for j=1:4
        X1(i,j)=X1(i,j)*(2^(4-j));%将每一行4个二进制数看成是一个4位的二进制数
    end
        source(i,1)=1+sum(X1(i,:));%%source(:,1)的每一行为X1矩阵中对应行的4位二进制数的值+1（值域1-16）
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
            3*d  -3*d];%mapping是一个16*2的矩阵，每一行对应着一个复数的实部和虚部
 for i=1:length(bitdata)/4
     qam_data(i,:)=mapping(source(i),:);%qam_data的每一行有两个元素，分别表示对应数据的QAM映射结果的实部和虚部
 end
 complex_qam_data=complex(qam_data(:,1),qam_data(:,2));%将qam_data每行的两个元素组合起来，得到复数列向量complex_qam_data，即为输出

end

