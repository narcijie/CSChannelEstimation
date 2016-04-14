function [Demodule_bit_symbol]=Demodu_QAM16(qam16_complex_symbols)
%===============================================
%16QAM的解调函数。由16QAM的复数信号恢复原二进制
%输入：qam16_complex_symbols 表示输入数据的列向量，每个元素为复数
%输出：Demodule_bit_symbol 表示输出结果的行向量，每个元素为二进制实数
%作者：杜捷    2015年4月6日
%===============================================
complex_symbols=reshape(qam16_complex_symbols,length(qam16_complex_symbols),1);
d=1;
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
            3*d  -3*d];%定义映射表
complex_mapping=complex(mapping(:,1),mapping(:,2));%complex_mapping为一个16*1的列向量，每个元素为一个复数
for i=1:length(qam16_complex_symbols);%对于每一个qam符号，将其与映射表中的16种情况一一比较
  for j=1:16;
      metrics(j)=abs(complex_symbols(i,1)-complex_mapping(j,1));%metrics行向量保存了16种情况的比较结果
  end
  [min_metric  decode_symbol(i)]= min(metrics) ;  %decode_symbol（i）为第i个QAM符号的16中映射中误差绝对值最小的映射的索引号  
end
 decode_bit_symbol=de2bi((decode_symbol-1)','left-msb');
 %decode_bit_symbole为一个n*4的矩阵，其每一行4个二进制数由decode_symbol转换而来
 Demodule_bit_symbol=reshape(decode_bit_symbol',1,length(qam16_complex_symbols)*4);
 %将decode_bit_symbole按从左到右，从上到下的顺序延展成行向量，然后输出

       

end

