%% 画图
h = h_ch(1 : ofdm.N);
H = fft(h);
%% 图1 冲激响应

figure(1)
figure_t = 0 : ofdm.N - 1;

stem(figure_t, abs(h_ch(1:ofdm.N)))
hold on
stem( sort(ESPRIT.Tau), 0.05 * ones( 1, MDL.L - 1),'r');

%% 图2 频率响应
figure(2)
stem(abs(H),'b.');
hold on
stem(abs(H_Sparse2(:,1)),'go');
stem(abs(H_Sparse(:,1)),'r.');

%% 图4  MSE
figure(4);
semilogy(chan.snrdBV,mean(loop.MSE_LS,1),'b-o')
hold on 
semilogy(chan.snrdBV,mean(loop.MSE_SP,1),'r.-')
hold on 
semilogy(chan.snrdBV,mean(loop.MSE_SP2,1),'g*-')
hold off
legend('LS估计','传统字典算法','改进字典算法');
grid on

%% 图5 BER
figure(5);
semilogy(chan.snrdBV,mean(loop.BER_LSE,1),'b-o')
hold on 
semilogy(chan.snrdBV,mean(loop.BER_SP,1),'r.-')
hold on 
semilogy(chan.snrdBV,mean(loop.BER_SP2,1),'g*-')
hold off
legend('LS估计','传统字典算法','改进字典算法');
grid on