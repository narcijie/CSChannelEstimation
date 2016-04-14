
   function [h]=cost_207(Sel_TypeOfChan,T,fd)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%             COST207 CHANNEL MODEL    %%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % h=channel handle output
   % T=sample period
   % fd=maximum doppler frequency
   % Sel_TypeOfChan=type of channel choice,there are four different types 
   %  specified in Cost207 Model including RA(Rural Area)(when 
   %  Sel_TypeOfChan=1),TU(Typical Urban)(when  Sel_TypeOfChan=2),
   %  BU(Bad Urban)(when Sel_TypeOfChan=3),HT(Hilly Terrian)(when
   %  Sel_TypeOfChan=4) all the types belong to multi-pathes channel
   % Project Name: Echo Canceller
% Author      : Chang Liu
% Index       : [Keywords of the project]
% Simulator   : Mathworks MATLAB 7.6(R14), SP3
% Date Created: 2008-12-19
% Last Update : 2008-12-19
% Version:    : 
% Description : [detailed functional description of the module/function]
% References  : [The citation list.]
% Dependency  : [The called functions or scripts.]
%
% Copyright 2008-2009, DSP EC R&D Group.
%
%% Revision History
% Rev. NO.    :
% Version     :
% Rev. Date   :
% Mender      :
% Description :
% References  : [The citation list.]
% Dependency  : [The called functions or scripts.]
%%
% Rev. NO.    :
% Version     :
% Rev. Date   :
% Mender      :
% Description :
% References  : [The citation list.]
% Dependency  : [The called functions or scripts.]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %construct tow types of Bigaussian doppler model in cost207 channel model
  TypeOfDop1_bigaussian = doppler.bigaussian('SigmaGaussian1', 0.05, ...
                          'SigmaGaussian2', 0.1, ...
                          'CenterFreqGaussian1', -0.8, ...
                          'CenterFreqGaussian2', 0.4, ...
                           'GainGaussian1', 1, ...
                           'GainGaussian2', 1/10);
  
  TypeOfDop2_bigaussian=doppler.bigaussian('SigmaGaussian1', 0.1, ...
                          'SigmaGaussian2', 0.15, ...
                          'CenterFreqGaussian1', 0.7, ...
                          'CenterFreqGaussian2', -0.4, ...
                           'GainGaussian1', 1, ...
                           'GainGaussian2', 1/10^1.5);  
                
    switch (Sel_TypeOfChan)
       
       case[1]      %choose RA TYPE,6 paths
         tau=[0 0.1 0.2 0.3 0.4 0.5].*1e-6;  %delay of the RA channel
         pdb=[0 -4 -8 -12 -16 -20];  %power of the RA channel
         %these parameters selected from the article "Performance of an echo
         %canceller and channel estimater for On-channel repeaters in
         %DVB_T/H networks" IEEE Transaction on broadcasting 2007.3
         k=0.91/0.41;  %K factor in ricianchan of RA
         h=ricianchan(T,fd,k,tau,pdb,0.7*fd);
            
        case[2]     %choose TU TYPE ,12 pathes
            tau=[0 0.2 0.4 0.6 0.8 1.2 ...
             1.4 1.8 2.4 3.0 3.2 5.0 ].*1e-6;  %delay of the TU channel
            pdb=[ -4 -3 0 -2 -3 -5 ... 
              -7 -5 -6 -9 -11 -10];   %power of the TU  channel
            h=rayleighchan(T,fd,tau,pdb);

          h.DopplerSpectrum=[doppler.jakes doppler.jakes doppler.jakes ...
                          TypeOfDop1_bigaussian TypeOfDop1_bigaussian...
                          TypeOfDop1_bigaussian TypeOfDop1_bigaussian ...
                          TypeOfDop1_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian];      
          % change the DopplerSpectrumType of h  
                    
      case[3]     %choose BU TYPE ,12 pathes
           tau=[0 0.2 0.4 0.8 1.6 2.2 ...
             3.2 5.0 6.0 7.2 8.2 10.0].*1e-6;  %delay of the TU channel
           pdb=[ -7 -3 -1 0 -2 -6 ... 
              -7 -1 -2 -7 -10 -15];   %power of the TU  channel
            h=rayleighchan(T,fd,tau,pdb);

          h.DopplerSpectrum=[doppler.jakes doppler.jakes doppler.jakes ...
                          TypeOfDop1_bigaussian TypeOfDop1_bigaussian...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian ...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian];      
           % change the DopplerSpectrumType of h  
                     
       case[4]     %choose HT TYPE ,12 pathes 
            tau=[0 0.2 0.4 0.8 1.6 2.0 ...
             2.4 15.0 15.2 15.8 17.2 20.0].*1e-6;  %delay of the TU channel
            pdb=[ -10 -8 -6 -4 0 0 ... 
              -4 -8 -9 -10 -12 -14];   %power of the TU  channel
            h=rayleighchan(T,fd,tau,pdb);

          h.DopplerSpectrum=[doppler.jakes doppler.jakes doppler.jakes ...
                          TypeOfDop1_bigaussian TypeOfDop1_bigaussian...
                          TypeOfDop1_bigaussian TypeOfDop2_bigaussian ...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian TypeOfDop2_bigaussian...
                          TypeOfDop2_bigaussian];      
           % change the DopplerSpectrumType of h  
       
      otherwise
           error('Unkown channel type!');
   end