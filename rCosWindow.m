function [ rcosw ] = rCosWindow( beta, Ts )
%===============================================
%定义升余弦窗函数
%输入：beta 滚降系数  ； Ts 包含循环前缀的OFDM符号的长度，为正偶数
%输出：rconw 表示输出结果的列向量
%作者：杜捷    2015年4月6日
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
rcosw = rcosw'; %变换为列向量
end

