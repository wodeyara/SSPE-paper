% Recreates the data for Figure 5 of the SSPE paper. 
%%
Fs = 1000;
timeAmt = 15;
L = Fs * timeAmt;
dt = 1/Fs;
ang_var2dev = @(v)(sqrt(-2*log(1-v)));

%%

% now running FIR + Hilbert
fNQ = Fs/2;
locutoff = 4;                               % Low freq passband = [4,8] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients



alpha = 1;
cnt2 = 1;
for SNR = [10]
%     initParams.sigmaFreqs = [.1/(SNR),.01*SNR];
    for cnt = 1:1000
        initParams = struct();
        initParams.freqs = [.1,6]; % drop first oscillator for the SNR = 10 case
        initParams.ampVec = [.99,.99];
        initParams.sigmaFreqs = [.5,.005];
        initParams.sigmaObs = .1;
        initParams.window = 5000;
        initParams.Fs = Fs;
        initParams.lowFreqBand = [4,8];
         %% generate data
        x1 = randn(L,1);
        xf1 = fft(x1);
        A = abs(xf1);
        phase = angle(xf1);

        df = 1.0 / (dt*length(x1));
        faxis = (0:(length(x1)/2))*df;
        faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
        oneOverf = 1.0 ./ faxis.^alpha;
        oneOverf(1)=1;

        Anew = sqrt((A.^2).*oneOverf');
        xf1new = Anew .* exp(1i*phase);
        x1new = (ifft(round(xf1new,8)))';
        x1new = x1new * 1000;

        x1 = randn(L,1);
        xf1 = fft(x1);
        A = abs(xf1);
        phase = angle(xf1);

        df = 1.0 / (dt*length(x1));
        faxis = (0:(length(x1)/2))*df;
        faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
        oneOverf = 1.0 ./ faxis.^alpha;
        oneOverf(1)=1;

        Anew = sqrt((A.^2).*oneOverf');
        xf1new = Anew .* exp(1i*phase);
        x2new = (ifft(round(xf1new,8)))';
        x2new = x2new * 1000;

        lowAct = filtfilt(filtwts,1,x2new);
        % SNR as std dev ratio of signal/noise
        SNR_start = std(lowAct)/std(x1new);
        data = (SNR/(SNR_start))*lowAct + x1new; % SNR = 1; sigmaFreqs [1,.05]
        data = data/std(data);
        truePhase = angle(hilbert((SNR/(SNR_start))*lowAct));
        
        simData(cnt2,cnt).data = data;
        simData(cnt2,cnt).truePhase = truePhase;
        
        %%

        [MK_phase,phaseBounds,fullX,phaseWidth,returnParams] = causalPhaseEM_MKmdl(data, initParams);
        MK_phase = reshape(MK_phase', size(MK_phase,1) * size(MK_phase,2),1);
        phaseWidth = reshape(phaseWidth', size(phaseWidth,1) * size(phaseWidth,2),1);
        phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
        
        if sum(returnParams.freqs>initParams.lowFreqBand(1) & returnParams.freqs<initParams.lowFreqBand(2))>1
            initParams.freqs = [6]; % drop first oscillator for the SNR = 10 case
            initParams.ampVec = [.99];
            initParams.sigmaFreqs = [.05];
            [MK_phase,phaseBounds,fullX,phaseWidth,returnParams] = causalPhaseEM_MKmdl(data, initParams);
            MK_phase = reshape(MK_phase', size(MK_phase,1) * size(MK_phase,2),1);
            phaseWidth = reshape(phaseWidth', size(phaseWidth,1) * size(phaseWidth,2),1);
            phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
        end
        
        simData(cnt2,cnt).SP_phase = MK_phase;
        simData(cnt2,cnt).SP_phaseBounds = phaseBounds;
        simData(cnt2,cnt).SP_phaseWidth = phaseWidth;
        simData(cnt2,cnt).SP_data = fullX;

        %%

        lowAct = filtfilt(filtwts,1,data);    
        lowAct_analytic = hilbert(lowAct);
        filt_phase = angle(lowAct_analytic);
        simData(cnt2,cnt).filterData = lowAct_analytic;
        simData(cnt2,cnt).filt_phase = angle(lowAct_analytic);      
        
        [phase_AW, ~, ~, analytic] = hilbert_transformer_causal(data', 1000,[4,8]);
        simData(cnt2,cnt).AW_phase = phase_AW;
        simData(cnt2,cnt).AW_data = analytic;

        %%

        %determine default parameters (median of optimized values)
        default_parameters = [];
        default_parameters.window_length = 750;%ceil(median(T.optim_window));
        default_parameters.filter_order = 192;%ceil(median(T.optim_filter_ord));
        default_parameters.edge = 64;% ceil(median(T.optim_edge));
        default_parameters.ar_order = 30;% ceil(median(T.optim_ar_ord));

        D = designfilt('bandpassfir', 'FilterOrder', default_parameters.filter_order, ...
            'CutoffFrequency1', 4, 'CutoffFrequency2', 8, 'SampleRate', Fs, 'DesignMethod', 'window');

        epochs = create_epochs_overlapping(data,Fs);
        epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);
        [zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
        simData(cnt2,cnt).zrenner_phase = zrenner_phase;
        simData(cnt2,cnt).estamp = estamp;

        circVarFIR(cnt2,cnt) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(filt_phase(5001:14e3) - truePhase(5001:14e3)))))));
        circVarAWMeth(cnt2,cnt) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(phase_AW(5001:14e3) - truePhase(5001:14e3)))))));
        circVarMKMeth(cnt2,cnt) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(MK_phase(5001:14e3)' - truePhase(5001:14e3)))))));
        circVarZrenner(cnt2,cnt) = rad2deg(ang_var2dev(1 - abs(mean(exp(1i*(zrenner_phase(4002:end) - truePhase(5001:14e3)))))));
    end
    cnt2 = cnt2 + 1;
end
%%

for SNR_i = 1:5
    for iter_i = 1:1000
        phaseWidths(SNR_i,iter_i) = mean(simData(SNR_i,iter_i).SP_phaseWidth(5001:14e3));
    end
end