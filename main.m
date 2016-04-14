
% clear
    clear                   % clear all variables
    close all               % close all figures
    clc                     % clear command window
%     randn('seed',1214)      % setting the seed for normal random numbers
%     rand('seed',12524)      % setting the seed for uniform random numbers
    
% ofdm parameters
    ofdm.N  =  128;                     % number of subcarriers   128
    ofdm.B  =  50;                       % number of block in each channel realization
    ofdm.M  =  4;                       % Modulation order (M=4)
    ofdm.T  =  5e-8;% 10ns              % OFDM sample time  1e-7
    ofdm.GI =  16;                      % length of gaurd interval  16
    ofdm.TS =  ofdm.N*ofdm.T;           % OFDM symbol time (not considering gaurd interval)
    ofdm.TT = (ofdm.N+ofdm.GI)*ofdm.T;  % OFDM symbol time (considering gaurd interval)
%     ofdm.PP =  1:10:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.Df =  5;
    ofdm.PP =  1:ofdm.Df:ofdm.N;              % Pilot locations in the subcarriers
    ofdm.DP =  setxor(1:ofdm.N,ofdm.PP);% Data  locations in the subcarriers
    ofdm.NP =  length(ofdm.PP);         % number subcarriers which carry data
    
% channel parameters
    chan.L      = 6;                                % number of channel taps
    chan.fd     = 0;%.1;                               % doppler in Hz
    chan.Nt     = 128;                              % Number of columns in the dictionary
    chan.Gain   = (0:1/(chan.Nt):1)*0;              % delay spread profile
    [~,chan.Delay]  = sort([0,rand(1,chan.Nt)]);    % generating random delay for each ta[
    chan.snrdB  = 15;                               % channel signal to noise ration
    chan.snrdBV = 10:2:30;                           % channel signal to noise ration for sweep
    
% loop parameters
%     loop.End1  = 1e2;                               % number of iterations
    loop.End1  = 1;                               % number of iterations
    loop.End2  = length(chan.snrdBV);               % length of inner loop
    loop.Sparse  = zeros(loop.End1,loop.End2);      % memory allocation for the BER using sparse method
    loop.BER_LSE     = zeros(loop.End1,loop.End2);      % memory allocation for the BER using LSE method
    loop.BER_SP     = zeros(loop.End1,loop.End2);
    loop.BER_SP2     = zeros(loop.End1,loop.End2);
    loop.MSE_SP  = zeros(loop.End1,loop.End2);
    loop.MSE_LS  = zeros(loop.End1,loop.End2);
    loop.MSE_SP2  = zeros(loop.End1,loop.End2);   
% building dictionary (please check different papers to learn how to build the dictionary)
%     chan.tau_p = linspace(0, ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./chan.Nt, chan.Nt);
% %     chan.tau_p = linspace(0, ofdm.TS, chan.Nt);
%     chan.Gamma = exp(-sqrt(-1)*2*pi.*repmat(  ((1:ofdm.N).') ,1,chan.Nt)    ./ ofdm.TS.* repmat(chan.tau_p,ofdm.N,1));
%     chan.tau_p2 = linspace(0, ofdm.GI*ofdm.T - ofdm.GI*ofdm.T./ 6, 6);
%     chan.Gamma2 = exp(-sqrt(-1)*2*pi.*repmat(  ((1:ofdm.N).') ,1,6)    ./ ofdm.TS.* repmat(chan.tau_p2,ofdm.N,1));
% fft matrix for LSE estimation (please check the following link for LSE channel estimation :
% http://www.mathworks.com/matlabcentral/fileexchange/13127-ofdm-lse-channel-estimation )

%     F = exp(2*pi*sqrt(-1)/ofdm.N .* meshgrid(0:ofdm.N-1,0:ofdm.N-1)...
%         .* repmat((0:ofdm.N-1)',[1,ofdm.N]));
    
% % Save the H
%     loop.SpH = zeros(loop.End1, loop.End2, ofdm.N);
%     loop.LsH = zeros(loop.End1, loop.End2, ofdm.N);
%     loop.Reh = zeros(loop.End1, loop.End2, ofdm.N);
%     
    loop.Chan_r = zeros(loop.End1, loop.End2, chan.L);
    loop.Chan_recog = zeros(loop.End1, loop.End2, chan.L);
%     loop.ReH = zeros(loop.End1, loop.End2, ofdm.N);
% Save MDL-L
    loop.MDLL = zeros(loop.End1, loop.End2);
%     loop.MDLy = zeros(12, loop.End2);
    hwait=waitbar(0,'正在计算');
    t_start = clock;
%% Loop 
for cnt1 = 1 :  loop.End1
    for cnt2 = 1 : loop.End2
        % loop parameters
        chan.snrdB = chan.snrdBV(cnt2);
        % Data generation
        data  = randi([0 ofdm.M-1],ofdm.N,ofdm.B);
        % modulation
        if ofdm.M == 4
            dataMod = qammod(data,ofdm.M)/sqrt(2);
        else
            error('Not defined')
        end
        
        % pilot insertion
        ofdm.Pilot = ones(ofdm.NP,ofdm.B);% or randsrc(ofdm.NP,ofdm.B,[-1 1]).*exp(-sqrt(-1)*pi*rand(ofdm.NP,ofdm.B));
        dataMod(ofdm.PP,:) = ofdm.Pilot;
        
        % ifft operation
        dataIFFT   = sqrt(ofdm.N)*ifft(dataMod);
        
        % adding gaurd interval
        dataIFFTGI = [dataIFFT((ofdm.N-ofdm.GI+1):ofdm.N,:);dataIFFT;];
        NDataGI = size(dataIFFTGI,1);
        
        % channel (rayleigh and gaussian noise)
%         ch = rayleighchan(ofdm.T, chan.fd, chan.tau_p(chan.Delay(1:chan.L)), chan.Gain(chan.Delay(1:chan.L)));
%         ch = cost_207(1, ofdm.T, chan.fd);      % 6 multipath
%         [ch loop.Chan_r(cnt1,cnt2,:)] = TestChannel(ofdm.T, 0, chan.L, floor(ofdm.N / ofdm.Df));
        [ch loop.Chan_r(cnt1,cnt2,:)] = TestChannel(ofdm.T, 0, chan.L, ofdm.GI);
%         dataChann = awgn(filter(ch,dataIFFTGI(:)),chan.snrdB );
        Pulse = [1 zeros(1, NDataGI - 1)]';
%         tData = dataIFFTGI(:);
        tData = [Pulse; dataIFFTGI(:)];
        dataChann1 = filter(ch,tData);
        h_ch = dataChann1(1 : NDataGI );
        dataChann2 = dataChann1(NDataGI + 1 : NDataGI * (ofdm.B + 1));
        dataChann = awgn(dataChann2, chan.snrdB);
        real_h = h_ch(1 : ofdm.N);
        real_H = fft(real_h);
        
%         realChannel = filter(ch,[1 zeros(1, ofdm.N - 1)]');


%         figure(3)
%         hold off
%         plot(abs(realChannel),'r');
%         loop.Reh(cnt1,cnt2,:) = filter(ch,[1 zeros(1, ofdm.N - 1)]');
%         loop.ReH(cnt1,cnt2,:) = fft(filter(ch,[1 zeros(1, ofdm.N - 1)]'),ofdm.N);
        % reshaping the signal
        dataChann = reshape(dataChann,ofdm.N+ofdm.GI,ofdm.B);
        
        % Guard interval removal
        dataRx     = dataChann((ofdm.GI+1):(ofdm.N+ofdm.GI),:);

        % ofdm demodulation
        dataRxFFT  =1/sqrt(ofdm.N)*fft(dataRx);
        
    
        %% LSE        
%         H_LSE = zeros(ofdm.N,ofdm.B);
%         for b = 1 : ofdm.B
% %              H_LSE(:,b) = ofdm.N/ofdm.NP * fft( ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b)),ofdm.N );
%             H_LSE(:,b) =  fft( ifft(dataRxFFT(ofdm.PP,b)./dataMod(ofdm.PP,b), ofdm.N),ofdm.N );
%              % 把一个只有ofdm.PP个长度的向量做ofdm.N长度的fft
%              % 相当于后面补0，再进行ofdm.N长度的fft
%             delta = abs(H_LSE(:,b) - real_H).^2;
%             loop.MSE_LS(cnt1,cnt2) =loop.MSE_LS(cnt1,cnt2) + mean(delta./ (abs(real_H).^2));
%         end
%         loop.MSE_LS(cnt1,cnt2) =loop.MSE_LS(cnt1,cnt2)./ ofdm.B;    
%         dataRxMod_LSE =  dataRxFFT(ofdm.DP,:)./H_LSE(ofdm.DP,:);             
%         dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);       
%         [~,BER_LSE] = biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M);
        
%         % saving the output
%         loop.Sparse(cnt1,cnt2) = BER_Sparse;
%         loop.LSE(cnt1,cnt2)    = BER_LSE;
%         
%         loop.SpH(cnt1,cnt2,:) = H_Sparse(:,1);
%         loop.LsH(cnt1,cnt2,:) = H_LSE(:,1);

        %%  LSE-传统
        H_LS = zeros(ofdm.N,ofdm.B);
        for b = 1 : ofdm.B
            H_LS(ofdm.PP,b) = dataRxFFT(ofdm.PP,b) ./ dataMod(ofdm.PP,b);             
            H_LS(ofdm.DP,b)=spline(ofdm.PP,H_LS(ofdm.PP,b),ofdm.DP);         
            delta = abs(H_LS(:,b) - real_H).^2;
            loop.MSE_LS(cnt1,cnt2) =loop.MSE_LS(cnt1,cnt2) + mean(delta./ (abs(real_H).^2));
        end
         loop.MSE_LS(cnt1,cnt2) =loop.MSE_LS(cnt1,cnt2)./ ofdm.B;  
        % 计算MSE
%         delta = abs(H_LS(:,b) - real_H);
%         loop.MSE_SP(cnt1,cnt2) = mean(delta./abs(real_H).^2);
        % 计算BER
        dataRxMod_LSE =  dataRxFFT(ofdm.DP,:) ./ H_LS(ofdm.DP,:);             
        dataRxDeMod_LSE = qamdemod(dataRxMod_LSE,ofdm.M);       
        loop.BER_LSE(cnt1,cnt2)= biterr(dataRxDeMod_LSE,data(ofdm.DP,:),ofdm.M) /  ( (ofdm.N - ofdm.NP) * ofdm.B );
        
        %% DFT方法
% %         H_DFT = zeros(ofdm.N,ofdm.B);
%         H_DFT = H_LSE;
%         h_DFT = ifft(H_DFT,ofdm.N) ./ sqrt(ofdm.N);
%         h_Sparse = ifft(H_Sparse,ofdm.N) ;
%         hold on
%         plot(abs(h_DFT),'g');
%         plot(abs(h_Sparse),'b');
%         hold off
%         for i = 20 : ofdm.N
%             h_DFT(i) = 0;
%         end
%         H_DFT = fft(h_DFT,ofdm.N);
%         dataRxMod_DFT =  dataRxFFT(ofdm.DP,:)./H_DFT(ofdm.DP,:);             
%         dataRxDeMod_DFT = qamdemod(dataRxMod_DFT,ofdm.M);       
%         [~,BER_DFT] = biterr(dataRxDeMod_DFT,data(ofdm.DP,:),ofdm.M);
%         loop.DFT(cnt1,cnt2) = BER_DFT;
%         

         %% 参数化信道估计 - MDL
        H_LS_P = zeros(ofdm.B, ofdm.NP);
        for b = 1 : ofdm.B
             H_LS_P(b, :) = dataRxFFT(ofdm.PP,b) ./ dataMod(ofdm.PP,b);
             % 导频处的LS估计
        end
        MDL.K = 5;
        MDL.I = ofdm.B;
        MDL.M = ofdm.NP;
        MDL.Q = zeros(MDL.M - MDL.K + 1, MDL.K - 1);
        MDL.R = zeros(MDL.M - MDL.K + 1, MDL.M - MDL.K + 1);
        J = rot90(eye(MDL.M - MDL.K + 1));
        for i = 1 : MDL.I
            for m = 1 : (MDL.M - MDL.K + 1)
                for n = 1 : (MDL.K)
                    MDL.Q(m,n) = H_LS_P(i, m + n - 1);
                end
            end
            MDL.R = MDL.R + (1 / (2 * MDL.K)) * (MDL.Q * MDL.Q' + J * conj(MDL.Q * MDL.Q') * J);
            % 求R矩阵的累加
        end
        MDL.R = MDL.R ./ MDL.I;        % 求相关阵R的平均值        
        MDL.lamda = eig(MDL.R);
        [MDL.V, MDL.D] = eig(MDL.R);   % MDL.D为特征向量矩阵，其每个列向量是一个特征向量 
        [MDL.lamda, MDL.lamdaIndex] = sort(MDL.lamda, 'descend');
        for t = 0 : (MDL.M - MDL.K)
            MDL.ms = 1;
            MDL.ps = 0;
            for tt = (t + 1) : (MDL.M - MDL.K + 1)
                MDL.ms = MDL.ms * MDL.lamda(tt);
                MDL.ps = MDL.ps + MDL.lamda(tt);
            end
            MDL.y(t + 1) = - MDL.I * (MDL.M - MDL.K + 1 - t) * ...
                log((MDL.ms ^ (1 / (MDL.M -MDL.K + 1 - t))) / (MDL.ps / (MDL.M - MDL.K + 1 - t))) + ...
                0.25 * t * (2 * (MDL.M - MDL.K + 1) - t + 1)* log(MDL.I);
        end
        [MDL.min, MDL.L] = min(MDL.y);
%         loop.MDLy(:,cnt2) = MDL.y;        
        loop.MDLL(cnt1, cnt2) = MDL.L - 1;
        
        %% ESPRIT medhod to acquire the initial multipath time delays
        MDL.L = 7;
        MDL.U = zeros(MDL.M - MDL.K + 1, MDL.M - MDL.K + 1);
        MDL.U(:) = MDL.V(:, MDL.lamdaIndex);       % (M - K + 1) * (M - K + 1)
        ESPRIT.U = MDL.U(:, 1 : (MDL.L - 1));      % (M - K + 1) * L
        ESPRIT.I1 = [eye(MDL.M - MDL.K) zeros(MDL.M - MDL.K, 1)];     % (M - K) * (M - K + 1)
        ESPRIT.I2 = [zeros(MDL.M - MDL.K, 1) eye(MDL.M - MDL.K)];     % (M - K) * (M - K + 1)
        ESPRIT.U1 = ESPRIT.I1 * ESPRIT.U;       % (M - K) * L
        ESPRIT.U2 = ESPRIT.I2 * ESPRIT.U;       % (M - K) * L
        
        ESPRIT.Phy = inv(ESPRIT.U1' * ESPRIT.U1) * ESPRIT.U1' * ESPRIT.U2;  % L * L
        ESPRIT.lamda = eig(ESPRIT.Phy);
        
        ESPRIT.Tau = zeros(1, MDL.L - 1);
        for i = 1 : MDL.L - 1
            ESPRIT.lamdaArg(i) = angle(ESPRIT.lamda(i)');    
            if ESPRIT.lamdaArg(i) < 0
                ESPRIT.lamdaArg(i) = 2 * pi + ESPRIT.lamdaArg(i);
            end            
            ESPRIT.Tau(i) = ESPRIT.lamdaArg(i) * ofdm.N / (2 * pi);
%             ESPRIT.Tau(i) = (angle(ESPRIT.lamda(i)') + pi) * ofdm.N / (2 * pi * ofdm.Df) * ofdm.T;
            ESPRIT.Tau(i) = mod(ESPRIT.Tau(i) + 1, 128) - 1;
            ESPRIT.Tau(i) = ESPRIT.Tau(i) ./ ofdm.Df;
        end
%         disp(ESPRIT.Tau);
        disp(sort(ESPRIT.Tau));
        loop.Chan_recog(cnt1, cnt2, :) = sort(ESPRIT.Tau);
            %% Sparse Channel estimation
        %动态字典算法---杜捷
        tau_rt = sort(ESPRIT.Tau);
        tau_len = length(tau_rt);
        nt_per_tau = floor(chan.Nt / tau_len);      % 每个多径的字典项数目
        tau_pp_1 = zeros(nt_per_tau, tau_len - 1);
        tau_lamda = 0.3;
        for i = 1 : tau_len - 1
            if tau_rt(i) < 0
                tau_rt(i) = 0;
            end
            if tau_rt(i) < 0.5
                tau_pp_1(:,i) = linspace( 0, (tau_rt(i) + tau_lamda) * ofdm.T, nt_per_tau);                
            else
                tau_pp_1(:,i) = linspace(  (tau_rt(i) - tau_lamda) * ofdm.T,  (tau_rt(i) + tau_lamda) * ofdm.T, nt_per_tau);
            end            
        end
        tau_pp_2 = linspace(  (tau_rt(tau_len) - tau_lamda) * ofdm.T,  (tau_rt(tau_len) + tau_lamda) * ofdm.T, chan.Nt - nt_per_tau * ( tau_len - 1));
        tau_pp_1 = reshape(tau_pp_1, 1, nt_per_tau * ( tau_len - 1) );
        tau_pp = [tau_pp_1 tau_pp_2]; 
        %动态字典算法---杜捷
        
        chan.tau_p = linspace(0,  ofdm.GI * ofdm.T - ofdm.GI*ofdm.T./chan.Nt ,  chan.Nt);       %ofdm.GI * ofdm.T就是最长可能的时间, 
        chan.tau_p2 = tau_pp;
        chan.Gamma = exp(-sqrt(-1)*2*pi.*repmat(  ((1:ofdm.N).') ,1,chan.Nt)    ./ ofdm.TS.* repmat(chan.tau_p,ofdm.N,1));
        chan.Gamma2 = exp(-sqrt(-1)*2*pi.*repmat(  ((1:ofdm.N).') ,1,chan.Nt)    ./ ofdm.TS.* repmat(chan.tau_p2,ofdm.N,1));
        H_Sparse = zeros(ofdm.N,ofdm.B);
        H_Sparse2 = zeros(ofdm.N,ofdm.B);
        lambda1 = ofdm.NP*10^(-chan.snrdB/10)/sum(abs(ch.pathGains));
        for b = 1 : ofdm.B
            y = dataRxFFT(ofdm.PP,b);
            A = chan.Gamma(ofdm.PP,:).*repmat(ofdm.Pilot(:,b),1,chan.Nt);           % chan.Nt: Number of columns in the dictionary
            A2 = chan.Gamma2(ofdm.PP,:).*repmat(ofdm.Pilot(:,b),1,chan.Nt);
%                 A = chan.Gamma2(ofdm.PP,:).*repmat(ofdm.Pilot(:,b),1,6);           % 修改了A的维度2015年12月3日
            cvx_begin quiet
                variable x(chan.Nt) complex
%                     variable x(6) complex
                    % sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
                    % 稀疏最优化方程（A矩阵根据字典矩阵G得出，y是接收到的数据）
                    % x是信道参数（即信道参数在G矩阵的列向量作为基时的坐标值）
                    minimize( quad_form(y-A*x, eye(ofdm.NP))+lambda1*norm(x,1) )
%                         minimize( quad_form(y-A*x, eye(ofdm.NP)) )
            cvx_end
            cvx_begin quiet
                variable x2(chan.Nt) complex
%                     variable x(6) complex
                    % sparse minimization formula (A is built from dictionary, y is received data and x is the channel coeff at pilot locations)
                    % 稀疏最优化方程（A矩阵根据字典矩阵G得出，y是接收到的数据）
                    % x是信道参数（即信道参数在G矩阵的列向量作为基时的坐标值）
                    minimize( quad_form(y-A2*x2, eye(ofdm.NP))+lambda1*norm(x2,1) )
%                         minimize( quad_form(y-A*x, eye(ofdm.NP)) )
            cvx_end
            % building channel at all location (simply from the dictionar)
            H_Sparse(:,b) = chan.Gamma*x;        
            H_Sparse2(:,b) = chan.Gamma2*x2; 
%         end
        % 计算MSE
        delta = abs(H_Sparse(:,b) - real_H).^2;
        loop.MSE_SP(cnt1,cnt2) =loop.MSE_SP(cnt1,cnt2) + mean(delta./ (abs(real_H).^2));
        delta = abs(H_Sparse2(:,b) - real_H).^2;
        loop.MSE_SP2(cnt1,cnt2) =loop.MSE_SP2(cnt1,cnt2) + mean(delta./ (abs(real_H).^2));
        end
        loop.MSE_SP(cnt1,cnt2) =loop.MSE_SP(cnt1,cnt2)./ ofdm.B;
        loop.MSE_SP2(cnt1,cnt2) =loop.MSE_SP2(cnt1,cnt2)./ ofdm.B;
        % 计算BER
        dataRxMod_Sparse =  dataRxFFT(ofdm.DP,:)./H_Sparse(ofdm.DP,:);             
        dataRxDeMod_Sparse = qamdemod(dataRxMod_Sparse,ofdm.M);       
        loop.BER_SP(cnt1,cnt2)= biterr(dataRxDeMod_Sparse,data(ofdm.DP,:),ofdm.M) /  ( (ofdm.N - ofdm.NP) * ofdm.B );
        
        dataRxMod_Sparse2 =  dataRxFFT(ofdm.DP,:)./H_Sparse2(ofdm.DP,:);             
        dataRxDeMod_Sparse2 = qamdemod(dataRxMod_Sparse2,ofdm.M);       
        loop.BER_SP2(cnt1,cnt2)= biterr(dataRxDeMod_Sparse2,data(ofdm.DP,:),ofdm.M) /  ( (ofdm.N - ofdm.NP) * ofdm.B );
        
        
        t_elapse = etime(clock,t_start);        
        t_total =   t_elapse / ( ( (cnt1 - 1)*loop.End2 + cnt2)/ (loop.End1* loop.End2) );
        t_remain = t_total - t_elapse;
        t_remain_min = floor(t_remain / 60);
        t_remain_sec = t_remain - t_remain_min*60;
        str=['已运行',num2str(( (cnt1 - 1)*loop.End2 + cnt2)/ (loop.End1* loop.End2)*100),'%已用时',num2str(floor(t_elapse)),'秒,剩余时间',num2str(t_remain_min) , '分' , num2str(floor(t_remain_sec)),'秒'];
        waitbar( ( (cnt1 - 1)*loop.End2 + cnt2)/ (loop.End1* loop.End2)  ,hwait,str);
    end
%     disp([num2str(round(cnt1/loop.End1*100)),'% has been done'])
end

     close(hwait); % 注意必须添加close函数
     
     t_spend = etime(clock,t_start);  
     t_spend_min = floor(t_spend / 60);
     t_spend_sec = floor(t_spend - t_spend_min*60);
     str = sprintf('计算共计用时%d分%d秒.', t_spend_min, t_spend_sec);
     disp(str);
%% Figure
plot_figure;
% 
% % disp(loop.MDLL);
% %     图1
% figure(1)
% figure_t = 0 : ofdm.N - 1;
% h = h_ch(1 : ofdm.N);
% H = fft(h);
% stem(figure_t, abs(h_ch(1:ofdm.N)))
% hold on
% stem( sort(ESPRIT.Tau), 0.05 * ones( 1, MDL.L - 1),'r');
% %%图2
% figure(2)
% stem(abs(H),'b.');
% hold on
% stem(abs(H_Sparse2(:,1)),'go');
% stem(abs(H_Sparse(:,1)),'r.');
% % stem(abs(H_LS(:,1)),'b*');
% 
% %%图3
% % figure(3)
% % stem(abs(h),'g.');
% % hold on
% % stem(abs( ifft(H_Sparse(:,1)) ),'ro');
% % stem(abs( ifft(H_LS(:,1)) ),'b*');
% %%图4 
% figure(4);
% semilogy(chan.snrdBV,mean(loop.MSE_LS,1),'b-o')
% hold on 
% semilogy(chan.snrdBV,mean(loop.MSE_SP,1),'r.-')
% hold on 
% semilogy(chan.snrdBV,mean(loop.MSE_SP2,1),'g*-')
% % % semilogy(chan.snrdBV,mean(loop.DFT,1),'g.-')
% % hold off
% % % legend('Sparse','LSE');%,'DFT')
% % grid on
% % 图5 BER
% figure(5);
% semilogy(chan.snrdBV,mean(loop.BER_LSE,1),'b-o')
% hold on 
% semilogy(chan.snrdBV,mean(loop.BER_SP,1),'r.-')
% hold on 
% semilogy(chan.snrdBV,mean(loop.BER_SP2,1),'g*-')
