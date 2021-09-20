% This simulation recreates the data for Figure 4 and Table 1 of the paper.
% To plot the figures look at plot_phaseEst_phaseReset.m
%%
clear 

Fs = 1000;
time = 10; 

indsToTestPhast = 1000:9000;
indsToTestPhase = setdiff(indsToTestPhast,1:2000);

initParams.freqs = [6];
initParams.Fs = 1000;
initParams.ampVec = [.99];
initParams.sigmaFreqs = [10];
initParams.sigmaObs = 1;
initParams.window = 2000;
initParams.lowFreqBand = [4,8];

fNQ = Fs/2;
locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

%determine default parameters (median of optimized values)
default_parameters = [];
default_parameters.window_length = 750;%ceil(median(T.optim_window));
default_parameters.filter_order = 192;%ceil(median(T.optim_filter_ord));
default_parameters.edge = 64;% ceil(median(T.optim_edge));
default_parameters.ar_order = 30;% ceil(median(T.optim_ar_ord));

D = designfilt('bandpassfir', 'FilterOrder', default_parameters.filter_order, ...
    'CutoffFrequency1', 4, 'CutoffFrequency2', 8, 'SampleRate', Fs, 'DesignMethod', 'window');

% timePointsSplice = [3500, 4750, 6500,8500];
timePointsSplice = [3500, 4750, 6500,8500];
indsToTestPhase = [];
for i = 1:4
    indsToTestPhase = [indsToTestPhase, timePointsSplice(i)+1:timePointsSplice(i)+167];
end
%%
ang_var2dev = @(v) sqrt(-2*log(1-v));

parfor iter = 1:1000
tic
          
    [pn] = make_pink_noise(1.5,time*Fs,1/Fs);
    % Make the data.
    V1 = (25).*cos(2*pi*(6).*[1/Fs:1/Fs:timePointsSplice(1)/Fs]);            
    V2 = (25).*cos(2*pi*(6).*[timePointsSplice(1)/Fs + 1/Fs:1/Fs:timePointsSplice(2)/Fs] + pi/2);
    V3 = (25).*cos(2*pi*(6).*[timePointsSplice(2)/Fs + 1/Fs:1/Fs:timePointsSplice(3)/Fs]);
    V4 = (25).*cos(2*pi*(6).*[timePointsSplice(3)/Fs + 1/Fs:1/Fs:timePointsSplice(4)/Fs] + pi/2);
    V5 = (25).*cos(2*pi*(6).*[timePointsSplice(4)/Fs + 1/Fs:1/Fs:time]);

    
    Vlo =  [V1,V2,V3,V4,V5];
     
    data = (Vlo) + (10)*pn; 
    truePhase = wrapTo2Pi([2*pi*(6).*[1/Fs:1/Fs:timePointsSplice(1)/Fs], ...
                        [2*pi*(6).*[timePointsSplice(1)/Fs + 1/Fs:1/Fs:timePointsSplice(2)/Fs] + pi/2], ...
                        [2*pi*(6).*[timePointsSplice(2)/Fs + 1/Fs:1/Fs:timePointsSplice(3)/Fs]], ...
                        [2*pi*(6).*[timePointsSplice(3)/Fs + 1/Fs:1/Fs:timePointsSplice(4)/Fs] + pi/2], ...
                        [2*pi*(6).*[timePointsSplice(4)/Fs + 1/Fs:1/Fs:time]]]);
    truePhase = truePhase';
    
    simData(iter).origData = data;
    simData(iter).truePhase = truePhase';
    
    epochs = create_epochs_overlapping(data,Fs);
%     [~, peak_SNR(iter)] = estimate_SNR(epochs, Fs, [4,8], []);
    % default parameters phastimate, fixed filter passband 8..13Hz
    epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);

    [phase,phaseBounds,fullX] = causalPhaseEM_MKmdl(data, initParams);
    SP_phase = reshape(phase', size(phase,1) * size(phase,2),1);
    phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
    fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));
    simData(iter).SP_phase = SP_phase;
    simData(iter).phaseBounds = phaseBounds;
    simData(iter).fullX = fullX;

    lowAct = filtfilt(filtwts,1,data);
    lowAct_analytic = hilbert(lowAct);   
    filt_phase = angle(lowAct_analytic);
    simData(iter).filt_phase = filt_phase;
    simData(iter).filtAnalytic = lowAct_analytic;
    
    [AW_phase,~,~, analytic] = hilbert_transformer_causal(data', 1000,[4,8]);
    simData(iter).AW_phase = AW_phase;
    simData(iter).AW_data = analytic;
      
    [zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
    simData(iter).zrenner_phase = zrenner_phase;
    simData(iter).zrenner_amp = estamp;
%     
%     err_SPcausal(iter) = abs(angle(mean(exp(1i*(truePhase(indsToTestPhase) - SP_phase(indsToTestPhase))))));
%     err_AW(iter) = abs(angle(mean(exp(1i*(truePhase(indsToTestPhase) - AW_phase(indsToTestPhase)')))));
%     err_FIR_Hilb(iter) = abs(angle(mean(exp(1i*(truePhase(indsToTestPhase) - filt_phase(indsToTestPhase)')))));
%     err_Zrenner(iter) = abs(angle(mean(exp(1i*(truePhase(indsToTestPhase) - zrenner_phase(1002:end)')))));
    
    var_SPcausal(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(indsToTestPhase) - SP_phase(indsToTestPhase)))))))
    var_AW(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(indsToTestPhase) - AW_phase(indsToTestPhase)'))))))
    var_FIR_Hilb(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(indsToTestPhase) - filt_phase(indsToTestPhase)'))))))
    var_Zrenner(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(indsToTestPhase) - zrenner_phase(indsToTestPhase-999)'))))))
    
%     var_SPcausal_2(iter) = 1 - abs(mean(exp(1i*(filt_phase(indsToTestPhase)' - SP_phase(indsToTestPhase)))));
%     var_AW_2(iter) = 1 - abs(mean(exp(1i*(filt_phase(indsToTestPhase)' - AW_phase(indsToTestPhase)'))));
%     var_FIR_Hilb_2(iter) = 1 - abs(mean(exp(1i*(filt_phase(indsToTestPhase)' - filt_phase(indsToTestPhase)'))));
%     var_Zrenner_2(iter) = 1 - abs(mean(exp(1i*(filt_phase(indsToTestPhase)' - zrenner_phase(1002:end)'))));
          
    
toc   
end
save('phaseReset_sines.mat')