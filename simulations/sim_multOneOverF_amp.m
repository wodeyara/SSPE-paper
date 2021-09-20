% simulates a 1/f spectrum for the amplitude, multiplies with sinusoid and
% has pink noise added. Very similar to Sines in Pink Noise simulation.

function sim_multOneOverF_amp()

% type = 'pink_noise';
time = 15;
Fs = 1000;
PAC = 0;
AAC = 0;

initParams.freqs = [.0001,6]; 
initParams.Fs = 1000;
initParams.ampVec = [.998,.999];
initParams.sigmaFreqs = [.001,.001];
initParams.sigmaObs = .01;
initParams.window = 5000;
initParams.lowFreqBand = [4,8];

% now running FIR + Hilbert
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

ang_var2dev = @(v)(sqrt(-2*log(1-v)));

for iter = 1:1000
    tic
   [data,truePhase] = generateData_R3sugg(Fs,time);
   simData(iter).origData = data;
   simData(iter).truePhase = truePhase;
   
   % first run causalEM
    try
    [SP_phase,phaseBounds,fullX,~,retParams] = causalPhaseEM_MKmdl_noSeg(data, initParams,1);
    
    simData(iter).SP_phase = SP_phase;
    simData(iter).SP_phaseBounds = phaseBounds;
    simData(iter).SP_data = fullX;
    simData(iter).SP_failed=0;
    err_SPcausal(iter) = mean(abs(SP_phase(2001:end) - truePhase(2001:end)))/(2*pi);
    catch
        simData(iter).SP_failed=1;
    end
    
    lowAct = filtfilt(filtwts,1,data);    
    lowAct_analytic = hilbert(lowAct);
    simData(iter).filterData = lowAct_analytic;
    simData(iter).filt_phase = angle(lowAct_analytic);

    err_FIR_Hilb(iter) = mean(abs(angle(lowAct_analytic(2001:end))' - truePhase(2001:end)))/(2*pi);
    
    % now running Anderson-Widge method
    [phase_AW,~,~, analytic] = hilbert_transformer_causal(data', 1000, [4,8]);
    simData(iter).AW_phase = phase_AW;
    simData(iter).AW_data = analytic;
    
    err_AW(iter) = mean(abs(phase_AW(2001:end)' - truePhase(2001:end)))/(2*pi);
    
    epochs = create_epochs_overlapping(data,Fs);
    epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);

    [zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
    simData(iter).zrenner_phase = zrenner_phase;
    simData(iter).zrenner_amp = estamp;
    
    err_zrenner(iter) = mean(abs(zrenner_phase(1002:end)' - truePhase(2001:14000)))/(2*pi);    
    
    var_SPcausal(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(5001:14000) - SP_phase(5001:14000)))))))
    var_AW(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(5001:14000) - phase_AW(5001:14000)'))))))
    var_FIR_Hilb(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(5001:14000) - angle(lowAct_analytic(5001:14000))'))))))
    var_Zrenner(iter) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(truePhase(5001:14000) - zrenner_phase(4002:end)'))))))

    toc       
end
   try 
    save(['causalPhaseEst_sim_MultoneOverF.mat'],'simData', 'err_AW', 'err_FIR_Hilb', 'err_SPcausal','err_zrenner',...
                                              'var_SPcausal','var_AW','var_FIR_Hilb','var_Zrenner');
   catch
       save('tmp.mat')
   end
end
    
function [data,truePhase] = generateData_R3sugg(Fs,timeAmt)  
    L = Fs*timeAmt; 
    dt = 1/Fs; 
    alpha = 4;
    freqOsc = 6;
    snr =.5; 

    x1 = randn(L,1);
    xf1 = fft(x1);

    df = 1.0 / (dt*length(x1));
    faxis = (0:(length(x1)/2))*df;
    faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
    oneOverf = 1.0 ./ (faxis.^(alpha/2));
    oneOverf(1)=0;

    Anew = (xf1.*oneOverf');
    x1new = real(ifft(round(Anew,8)));
    x1new = (x1new/std(x1new));        
    x1new = x1new + abs(min(x1new));

    [pn] = make_pink_noise(1.5,L,1/Fs);
    pn = pn/std(pn);
    lowAct = x1new'.*cos(2*pi*freqOsc*[1/Fs:1/Fs:timeAmt]);
    SNRstart = std(lowAct)/std(pn);
    data =  (snr/SNRstart)*lowAct + pn;
    data = data/std(data);
    truePhase = wrapToPi(2*pi*freqOsc*[1/Fs:1/Fs:timeAmt])';
%     truePhase = angle(hilbert(20*x1new'.*cos(2*pi*freqOsc*[1/Fs:1/Fs:timeAmt])))';
end
%%

% function [data,truePhase] = generateData(Fs,timeAmt)  
% 
% L = Fs*timeAmt;
% dt = 1/Fs; 
% alpha =1;
% 
% x1 = randn(L,1);
% xf1 = fft(x1);
% A = abs(xf1);
% phase = angle(xf1);
% 
% df = 1.0 / (dt*length(x1));
% faxis = (0:(length(x1)/2))*df;
% faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
% oneOverf = 1.0 ./ faxis.^alpha;
% oneOverf(1)=1;
% 
% Anew = sqrt((A.^2).*oneOverf');
% xf1new = Anew .* exp(1i*phase);
% x1new = (ifft(round(xf1new,8)))';
% x1new = (x1new/std(x1new));        
%       
% x1 = randn(L,1);
% xf1 = fft(x1);
% A = abs(xf1);
% phase = angle(xf1);
% 
% df = 1.0 / (dt*length(x1));
% faxis = (0:(length(x1)/2))*df;
% faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
% oneOverf = 1.0 ./ faxis.^alpha;
% oneOverf(1)=1;
% 
% phase(1) = 0;
% Anew = sqrt(A.^2.*oneOverf');
% Anew(1) = .5e4;
% xf1new = Anew .* exp(1i*phase);
% x2new = (ifft(round(xf1new,8)))';
% x2new = abs((x2new/std(x2new)))+.1;
% data = x2new.*cos(2*pi*6*[1/Fs:1/Fs:timeAmt]) + 2*x1new;
% truePhase = angle(hilbert(cos(2*pi*(6).*[1/Fs:1/Fs:timeAmt])))';
%     
% end