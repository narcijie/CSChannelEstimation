function [Demodule_bit_symbol]=Demodu_QAM16(qam16_complex_symbols)
%===============================================
%16QAM�Ľ����������16QAM�ĸ����źŻָ�ԭ������
%���룺qam16_complex_symbols ��ʾ�������ݵ���������ÿ��Ԫ��Ϊ����
%�����Demodule_bit_symbol ��ʾ����������������ÿ��Ԫ��Ϊ������ʵ��
%���ߣ��Ž�    2015��4��6��
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
            3*d  -3*d];%����ӳ���
complex_mapping=complex(mapping(:,1),mapping(:,2));%complex_mappingΪһ��16*1����������ÿ��Ԫ��Ϊһ������
for i=1:length(qam16_complex_symbols);%����ÿһ��qam���ţ�������ӳ����е�16�����һһ�Ƚ�
  for j=1:16;
      metrics(j)=abs(complex_symbols(i,1)-complex_mapping(j,1));%metrics������������16������ıȽϽ��
  end
  [min_metric  decode_symbol(i)]= min(metrics) ;  %decode_symbol��i��Ϊ��i��QAM���ŵ�16��ӳ����������ֵ��С��ӳ���������  
end
 decode_bit_symbol=de2bi((decode_symbol-1)','left-msb');
 %decode_bit_symboleΪһ��n*4�ľ�����ÿһ��4������������decode_symbolת������
 Demodule_bit_symbol=reshape(decode_bit_symbol',1,length(qam16_complex_symbols)*4);
 %��decode_bit_symbole�������ң����ϵ��µ�˳����չ����������Ȼ�����

       

end

