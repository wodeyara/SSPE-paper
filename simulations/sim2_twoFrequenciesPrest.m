% Recreates the data for Figure 

%%
% sim2
% Here we show that when there are two oscillations in the low freqeucny
% band we struggle to estimate the right phase, however the SP/MK approach
% does great at pulling this out.


% lets test both the hilbert transformer and the current method jointly:
% forcing data sent to both to be the same.

clear 

Fs  = 1000;
time= 15;
freq = 6;
cnt = 1;
mult = [0:4,5,4:-1:0];

ang_var2dev = @(v)(sqrt(-2*log(1-v)));

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

timeDelayAWMeth = zeros (10,11);
timeDelayMKMeth = zeros(10,11);
timeDelayFIR = zeros(10,11);
timeDelayZrenner = zeros(10,11);

circVarFIR = zeros(10,11);
circVarAWMeth = zeros(10,11);
circVarMKMeth = zeros(10,11);
circVarZrenner = zeros(10,11);

for secOscDist = [-5:5]
    tic
    ampOfSecOsc = [0.2:.2:2];
    parfor cnt2 = 1:length(ampOfSecOsc)      
        initParams = struct();
        data = (25).*cos(2*pi*(freq)*[1/Fs:1/Fs:time]) +...
            ampOfSecOsc(cnt2)*(25).*cos(2*pi*(freq+secOscDist)*[1/Fs:1/Fs:time] + pi/3) ...
            + .5*randn(1,time*Fs);
        
        truePhase =wrapToPi(2*pi*(freq)*[1/Fs:1/Fs:time]);
        
        lowAct = filtfilt(filtwts,1,data);
        lowAct_analytic = hilbert(lowAct);   
        filt_phase = angle(lowAct_analytic);

%         [init_f, init_a,init_sigma,R0] = initializeParams(data,2,1,5000);

        initParams.freqs = [5,6];
        initParams.Fs = 1000;
        initParams.ampVec = [.99,.99];
        initParams.sigmaFreqs = [10,20];
        initParams.sigmaObs = 1;
        initParams.window = 5000;
        initParams.lowFreqBand = [5.25,6.75];

        [phase,phaseBounds, fullX] = causalPhaseEM_MKmdl(data, initParams);

        phase = reshape(phase', size(phase,1) * size(phase,2),1);
        phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
        fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));
% 
        [phase_AW, ~, ~, analytic] = hilbert_transformer_causal(data', 1000,[4,8]);
        
        epochs = create_epochs_overlapping(data,Fs);
        epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);
        [zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
        
        timeDelayFIR(cnt2,cnt) = angle(mean(exp(1i*(filt_phase(5001:14e3) - truePhase(5001:14e3)))))*(180/pi);
        timeDelayAWMeth(cnt2,cnt) = angle(mean(exp(1i*(phase_AW(5001:14e3) - truePhase(5001:14e3)))))*(180/pi);
        timeDelayMKMeth(cnt2,cnt) = angle(mean(exp(1i*(phase(5001:14e3)' - truePhase(5001:14e3)))))*(180/pi);
        timeDelayZrenner(cnt2,cnt) = angle(mean(exp(1i*(zrenner_phase(4002:end) - truePhase(5001:14e3)))))*(180/pi);

        % circular variance
        circVarFIR(cnt2,cnt) = 1 - abs(mean(exp(1i*(filt_phase(5001:14e3) - truePhase(5001:14e3)))));
        circVarAWMeth(cnt2,cnt) = 1 - abs(mean(exp(1i*(phase_AW(5001:14e3) - truePhase(5001:14e3)))));
        circVarMKMeth(cnt2,cnt) = 1 - abs(mean(exp(1i*(phase(5001:14e3)' - truePhase(5001:14e3)))));
        circVarZrenner(cnt2,cnt) = 1 - abs(mean(exp(1i*(zrenner_phase(4002:end) - truePhase(5001:14e3)))));
    end
    toc
    cnt = cnt + 1;
end

circVarFIR = rad2deg(ang_var2dev(circVarFIR));
circVarAWMeth = rad2deg(ang_var2dev(circVarAWMeth));
circVarMKMeth = rad2deg(ang_var2dev(circVarMKMeth));
circVarZrenner = rad2deg(ang_var2dev(circVarZrenner));
%%
circVarFIR(:,6) = NaN;
circVarMKMeth(:,6) = NaN;
circVarZrenner(:,6) = NaN;
circVarAWMeth(:,6) = NaN;

subplot(221)
imagesc([0.2:.2:2], 6+[-5:5],circVarFIR')
caxis([0,120])
% set(gca,'YTick', 1:10, 'YTickLabel', 6+[-5:2:-1,1:2:5])
set(gca,'Fontsize', 14)
xlabel('Proportion of Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
tmp = gray;
colormap(flipud(tmp))
c = colorbar;
c.Label.String = 'Circular Standard Deviation';
 title('acausal FIR', 'Fontsize', 16)

subplot(222)
imagesc([0.2:.2:2], 6+[-5:5],circVarAWMeth')
caxis([0,120])
% set(gca,'YTick', 1:2:11, 'YTickLabel', 6+[-5:2:-1,1:2:5])
set(gca,'Fontsize', 14)
xlabel('Proportion of Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
% colormap(gray)
c = colorbar;
c.Label.String = 'Circular Standard Deviation';
 title('Blackwood et al. 2018', 'Fontsize', 16)


subplot(223)
imagesc([0.2:.2:2], 6+[-5:-1,1:5],circVarZrenner')
caxis([0,120])
% set(gca,'YTick', 1:2:11, 'YTickLabel', 6+[-5:2:-1,1:2:5])
set(gca,'Fontsize', 14)
xlabel('Proportion of Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
% colormap(gray)
c = colorbar;
c.Label.String = 'Circular Standard Deviation';
 title('Zrenner et al. 2020', 'Fontsize', 16)

subplot(224)

imagesc([0.2:.2:2], 6+[-5:5],circVarMKMeth')
caxis([0,120])
% set(gca,'YTick', [1:2:9,10], 'YTickLabel', 6+[-5:2:-1,1:2:5])
set(gca,'Fontsize', 14)
xlabel('Proportion of Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
% colormap(gray)
c = colorbar;
c.Label.String = 'Circular Standard Deviation';
 title('SSPE', 'Fontsize', 16)

%% testing 
% 
% Fs  = 1000;
% time= 10;
% freq = 6;
% cnt = 1;
% pinkNoise = dsp.ColoredNoise('pink',Fs,10);
% 
% ampOfSecOsc = 10;
% secOscFreq = 0.5;
% 
% data = 10*cos(2*pi*(freq)*[1/Fs:1/Fs:time]) + 10*ampOfSecOsc*cos(2*pi*(secOscFreq)*[1/Fs:1/Fs:time] + pi/4) ...
%             + reshape(step(pinkNoise),1,Fs*time);
% realPhase = wrapToPi((2*pi*(freq)*[1/Fs:1/Fs:time]));
% 
% params.Fs = 1000;
% freqsTested =  0.5:0.5:15;
% freqLoc = 1;
% for i = freqsTested
%     [~,Fstat(freqLoc)] = mtpowerandfstatc(data,params,i);
%     freqLoc = freqLoc + 1;
% end
% K = 5; % this is the default in mtpowerandfstatc
% pval = (1/length(freqsTested)) * 0.01; %bonferroni corrected and single sided test - no?
% freqStart = freqsTested(Fstat>finv(1 - pval,2,2*K-2));
% 
% [omega,newAllX,newAllP] = fitSoulatModel_multSines(data,freqStart,1000);
% [~,ind] = min(abs(omega-6));
% estPhaseSoulatModel = angle(newAllX((ind-1)*2 +1,:) + 1i*newAllX((ind-1)*2 +2,:));
% 
% [phase, estimate_mask, analytic] = hilbert_transformer_phase(data, 500);
% % getting metrics of interest
% dsRealPhase = wrapToPi((2*pi*(freq)*[1/1000:1/1000:time])); %1000 is downsampled sampling frequency
% 
% (mean((abs(estPhaseSoulatModel' - dsRealPhase'))/(2*pi) * (1/freq) *1000))


%%
% old approach

%         tic
%         params.Fs = 1000;
%         freqsTested =  1:200;
%         freqLoc = 1;
%         for i = freqsTested
%             [~,Fstat(freqLoc)] = mtpowerandfstatc(data,params,i);
%             freqLoc = freqLoc + 1;
%         end
%         K = 5; % this is the default in mtpowerandfstatc
%         pval = (1/length(freqsTested)) * 0.01; %bonferroni corrected and single sided test - no?
%         freqStart = freqsTested(Fstat>finv(1 - pval,2,2*K-2));
%         toc
%         
%          [~,ind] = min(abs(omega-6));
%         estPhaseSoulatModel = angle(newAllX((ind-1)*2 +1,:) + 1i*newAllX((ind-1)*2 +2,:));
