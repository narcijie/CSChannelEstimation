% TR-MIMO-Sparse Channel
%% Clean up
clear,clc
close all
%% Tx & Rx number, lx length
Nt = 2;% Tx
Nr = 10;% Rx
lx_pilot = 150;% X length
lx_unknown = 1e3 - lx_pilot;
lx = lx_pilot + lx_unknown;
%% Constellation
map = 1/sqrt(2) * [1+1i;-1+1i;-1-1i;1-1i;];
%% Sparse channel
L = 29;
lh = L + 1;
valNum = 5;
temp = randperm(lh);
chanIndex = temp(1:valNum);
hRaw = exp((valNum:-1:1)/valNum);
hProfile = zeros(1,lh);
for iAssignValue = 1:valNum
    hIndex = chanIndex(iAssignValue);
    hProfile(hIndex) = hRaw(iAssignValue);
end
hProfile_n = hProfile./norm(hProfile);
%% SNR
SNR_dB = 0:1:20;
nTrial = (1:length(SNR_dB)) * 1e1;
BER_aver = zeros(length(SNR_dB),Nt);
BER_aver_Ls = zeros(length(SNR_dB),Nt);
BER_aver_Omp = zeros(length(SNR_dB),Nt);
%% Iteration
for iSNR = 1:length(SNR_dB)
    iSNR,
    BER_sum = 0;
    BER_sum_Ls = 0;
    BER_sum_Omp = 0;
    for iTrial = 1:nTrial(iSNR)
        %% Signal and noise generate, transmit through AWGN channel
        tx = zeros(Nt,lx);
        hNorm = zeros(lh,Nr,Nt);
        noise = zeros(Nr,lh+lx-1);
        ySingle = zeros(lh+lx-1,Nr,Nt);
        yNoNoise = zeros(lh+lx-1,Nr);
        % Generate tx and channel
        for iTx = 1:Nt
            tx(iTx,:) = 1 / sqrt(2) * ((randi(2,1,lx) - 1.5) * 2 + (randi(2,1,lx) - 1.5) * 2 * 1i);
            for iRx = 1:Nr
                % Get channel, noise, and conv noise with tx
                hNorm(:,iRx,iTx) = 1/sqrt(2)*hProfile_n.'.*(randn(lh,1)+randn(lh,1)*1i);% randomly generate h
                ySingle(:,iRx,iTx) = conv(hNorm(:,iRx,iTx),tx(iTx,:));
            end
            yNoNoise = yNoNoise + ySingle(:,:,iTx);
        end
        % Generate noise
        n_Null = 10 ^ (-SNR_dB(iSNR) / 10);
        for iRx = 1:Nr
            noise(iRx,:) = 1 / sqrt(2) * sqrt(n_Null) * (randn(1,lh+lx-1)+randn(1,lh+lx-1)*1i);
        end
        z = noise.';
        % Get output
        y = yNoNoise + z;
        %% ChanEST-LS
        yLs = zeros(lh+lx-1,Nr,Nt);
        hHatLs = zeros(lh,Nr,Nt);
        xMatrix = zeros(lx_pilot,lh);
        for iTx = 1:Nt
            for iColumn = lh:-1:1
                xMatrix(:,iColumn) = tx(iTx,lh-iColumn+1:lh-iColumn+lx_pilot).';
            end
            yLs(:,:,iTx) = ySingle(:,:,iTx) + z;
            hHatLs(:,:,iTx) = (xMatrix' * xMatrix) \ xMatrix' * yLs(lh:lx_pilot+lh-1,:,iTx);
        end
        %% ChanEST-OMP
        magic = 1.2;
        hHatOmp = zeros(lh,Nr,Nt);
        yOmp = zeros(lh+lx-1,Nr,Nt);
        sparcity = ceil(valNum * magic);
        for iTx = 1:Nt
            for iColumn = lh:-1:1
                xMatrix(:,iColumn) = tx(iTx,lh-iColumn+1:lh-iColumn+lx_pilot).';
            end
            yOmp(:,:,iTx) = ySingle(:,:,iTx) + z;
            indx = zeros(1,sparcity);
            for iRx = 1:Nr
                residual = yOmp(lh:lx_pilot+lh-1,iRx,iTx);
                for iter = 1:sparcity
                    inProj = zeros(1,lh);
                    for iColumn = 1:lh
                        inProj(iColumn) = dot(residual,xMatrix(:,iColumn)) / norm(xMatrix(:,iColumn))^2;
                    end
                    [maxVal,pos] = max(abs(inProj));
                    indx(iter) = pos;
                    Csp = xMatrix(:,indx(1:iter));
                    hHat = (Csp'* Csp) \ Csp'* yOmp(lh:lx_pilot+lh-1,iRx,iTx);
                    residual = yOmp(lh:lx_pilot+lh-1,iRx,iTx) - Csp * hHat;
                end
                for ivalue = 1:length(indx)
                    hHatOmp(indx(ivalue),iRx,iTx) = hHat(ivalue);
                end
            end
        end
        %% Decode by time-revrsal
        [BER, xDecode] = decodeTR(1,Nt,Nr,lx,lx_unknown,L,hNorm,y,map,tx,lh,lx_pilot);
        BER_sum = BER_sum + BER;
        [BERLs, xDecodeLs] = decodeTR(2,Nt,Nr,lx,lx_unknown,L,hHatLs,y,map,tx,lh,lx_pilot);
        BER_sum_Ls = BER_sum_Ls + BERLs;
        [BEROmp, xDecodeOmp] = decodeTR(3,Nt,Nr,lx,lx_unknown,L,hHatOmp,y,map,tx,lh,lx_pilot);
        BER_sum_Omp = BER_sum_Omp + BEROmp;
%         [BERLs2, xDecodeLs2] = decodeTRchanEST(2,Nt,Nr,lx,lx_unknown,L,hHatLs,y,map,tx,lh,lx_pilot);
%         BER_sum_Ls = BER_sum_Ls + BERLs;
%         [BEROmp2, xDecodeOmp2] = decodeTRchanEST(3,Nt,Nr,lx,lx_unknown,L,hHatOmp,y,map,tx,lh,lx_pilot);
%         BER_sum_Omp = BER_sum_Omp + BEROmp;
    end
    BER_aver(iSNR,:) = (BER_sum / nTrial(iSNR)).';
    BER_aver_Ls(iSNR,:) = (BER_sum_Ls / nTrial(iSNR)).';
    BER_aver_Omp(iSNR,:) = (BER_sum_Omp / nTrial(iSNR)).';
end
%% Plot
semilogy(SNR_dB,BER_aver(:,1),'--ro');
hold on
semilogy(SNR_dB,BER_aver(:,2),'--bo');
hold on
semilogy(SNR_dB,BER_aver_Ls(:,1),'--r+');
hold on
semilogy(SNR_dB,BER_aver_Ls(:,2),'--b+');
hold on
semilogy(SNR_dB,BER_aver_Omp(:,1),'--rx');
hold on
semilogy(SNR_dB,BER_aver_Omp(:,2),'--bx');
grid on
title('Bit Error Rate for different Tr1 and Tr2');
xlabel('SNR(dB)');
ylabel('Bit Error Rate');
legend('Tr1','Tr2','Tr1-LS EST','Tr2-LS EST','Tr1-OMP EST','Tr2-OMP EST','Location','NorthEast');
