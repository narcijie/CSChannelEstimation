function y=Single_Path_Flat_Rayleigh_Fading(fm,fs,L)  
%**************************************************************************
% 此函数功能——输出单径平坦瑞利衰落的复高斯序列,且产生的该复高斯序列具有经
% 典的多普勒功率谱。其中fm表示最大多普勒频移,fs表示采样频率,即该复高斯序列是
% 以采样频率fs对一连续的复高斯随机过程采样得到的,L表示该复高斯序列的长度。此
% 外输出的该复高斯序列的均值为0,包络的均方值(包络的平均功率)为1;且该复高斯序               
% 列的实部与虚部独立不相关,且正交。
% 此复高斯序列的产生原理是基于典型的Jakes模型,即两路独立的高斯白噪声分别通过
% 一个窄带的多普勒滤波器。
%**************************************************************************

N1=2^12;
N2=2^18;
T0=(N1-0.5)/fm;
fs0=N2/T0;
K1=round(fs0/fs);
K2=round(fs/fs0);
if fs0>fs
    L0=L*K1;
else
    L0=ceil(L/K2);
end
n=0:(N1-1);
h1=1/sqrt(2*pi*fm)./(1-(n/(N1-0.5)).^2).^(1/4);
h2=h1(2:N1);
h3=fliplr(h2);
h4=zeros(1,N2-2*N1+1);
h5=[h1 h4 h3];
h6=real(ifft(h5));
h7=h6(2:N2/2);
h8=h6(1:N2/2);
h9=fliplr(h7);
h=[h9 h8];
x1=randn(1,N2-1+L0);
x2=randn(1,N2-1+L0);
L1=2*N2-3+L0;
H=fft(h,L1);
X1=fft(x1,L1);
X2=fft(x2,L1);
y1=sqrt(fs0)*real(ifft(H.*X1));
y2=sqrt(fs0)*real(ifft(H.*X2));
if fs0>fs
    y3=y1(N2-1:N2+L0-2);
    y4=y2(N2-1:N2+L0-2);
    y5=downsample(y3,K1);
    y6=downsample(y4,K1);
else
    y3=y1(N2-1:N2+L0-1);
    y4=y2(N2-1:N2+L0-1);
    x=1:K2:(L0+1)*K2;
    xi=1:L0*K2+1;
    y5=interp1(x,y3,xi,'linear');
    y6=interp1(x,y4,xi,'linear');
end
y7=y5(1:L);
y8=y6(1:L);
y=y7+j*y8;





