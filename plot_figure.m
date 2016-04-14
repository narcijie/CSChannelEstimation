%% ��ͼ
h = h_ch(1 : ofdm.N);
H = fft(h);
%% ͼ1 �弤��Ӧ

figure(1)
figure_t = 0 : ofdm.N - 1;

stem(figure_t, abs(h_ch(1:ofdm.N)))
hold on
stem( sort(ESPRIT.Tau), 0.05 * ones( 1, MDL.L - 1),'r');

%% ͼ2 Ƶ����Ӧ
figure(2)
stem(abs(H),'b.');
hold on
stem(abs(H_Sparse2(:,1)),'go');
stem(abs(H_Sparse(:,1)),'r.');

%% ͼ4  MSE
figure(4);
semilogy(chan.snrdBV,mean(loop.MSE_LS,1),'b-o')
hold on 
semilogy(chan.snrdBV,mean(loop.MSE_SP,1),'r.-')
hold on 
semilogy(chan.snrdBV,mean(loop.MSE_SP2,1),'g*-')
hold off
legend('LS����','��ͳ�ֵ��㷨','�Ľ��ֵ��㷨');
grid on

%% ͼ5 BER
figure(5);
semilogy(chan.snrdBV,mean(loop.BER_LSE,1),'b-o')
hold on 
semilogy(chan.snrdBV,mean(loop.BER_SP,1),'r.-')
hold on 
semilogy(chan.snrdBV,mean(loop.BER_SP2,1),'g*-')
hold off
legend('LS����','��ͳ�ֵ��㷨','�Ľ��ֵ��㷨');
grid on