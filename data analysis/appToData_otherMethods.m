% Given a vector of data and an intent to estimate phase for a theta band
% this code uses alternative real-time estimation techniques for the phase
% estimation.
%%
lowFreq = 4;
highFreq = 11;

% now running FIR + Hilbert
fNQ = Fs/2;
locutoff = lowFreq;                               % Low freq passband = [4,7] Hz.
hicutoff = highFreq;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

lowAct = filtfilt(filtwts,1,data);    
lowAct_analytic = hilbert(lowAct);
filterData = lowAct_analytic;
filt_phase = angle(lowAct_analytic);

%determine default parameters (median of optimized values)
default_parameters = [];
default_parameters.window_length = 750;%ceil(median(T.optim_window));
default_parameters.filter_order = 192;%ceil(median(T.optim_filter_ord));
default_parameters.edge = 64;% ceil(median(T.optim_edge));
default_parameters.ar_order = 30;% ceil(median(T.optim_ar_ord));

D = designfilt('bandpassfir', 'FilterOrder', default_parameters.filter_order, ...
    'CutoffFrequency1', lowFreq, 'CutoffFrequency2', highFreq, 'SampleRate', Fs, 'DesignMethod', 'window');
    
epochs = create_epochs_overlapping(data,Fs);
epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);

[zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
zrenner_phase = [zeros(1,999),zrenner_phase,zeros(1,1000)];

[phase_AW,~,~, analytic] = hilbert_transformer_causal(data', 1000, [lowFreq,highFreq]);

%% estimate phase and phase Width
initParams = allInitParams{41};
[phase, ~,~,phaseWidth] = runSSPE_KalmanFilter(data(50*Fs + 1:end),initParams);

%%
cnt = 1; truePhase = [];
for i = 46:length(allStateVec)
    trueX = allStateVec{i};
    tmpParamTrue = allInitParams{i};
    lowFreqLoc = find(tmpParamTrue.freqs>4 & tmpParamTrue.freqs<11);
    if ~isempty(lowFreqLoc) & length(lowFreqLoc) ==1
        truePhase(cnt,:) =  angle(trueX(lowFreqLoc*2 - 1,5001:6000) + trueX(lowFreqLoc*2,5001:6000)*1i);   
    end
    cnt = cnt +1;
end
 truePhase = reshape(truePhase', size(truePhase,1) * size(truePhase,2),1);
 %%

redData = data(50*Fs + 1:length(truePhase) + 50*Fs);
% phaseEst(:,1) = phase_AW(50*Fs + 1:length(truePhase) + 50*Fs);
% phaseEst(:,2) = filt_phase(50*Fs + 1:length(truePhase) + 50*Fs);
% phaseEst(:,3) = zrenner_phase(50*Fs + 1:length(truePhase) + 50*Fs);
% phaseEst(:,4) = phaseCredLim(1:length(truePhase));

for i = 1:length(truePhase)/Fs
    tmpDataOneSec = redData((i-1) * Fs + 1:i*Fs);
    tmpTruePhase = truePhase((i-1) * Fs + 1:i*Fs);
    tmpPhaseEst = phaseEst((i-1) * Fs + 1:i*Fs,:);
    
    if sum(tmpTruePhase) == 0
        continue
    end
    for j = 1:3
        trueError(i,j) = rad2deg(ang_var2dev(1-abs(mean(exp(1i*(tmpPhaseEst(:,j) - tmpTruePhase))))));
    end
end

trueError(trueError==0) = NaN;
trueErrorCredLim(trueErrorCredLim==0)  = NaN;
trueError(:,4) = trueErrorCredLim;
trueError = trueError(:,[2,1,3,4]);
%%
% general graphics, this will apply to any figure you open
% (groot is the default figure object).
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 8, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 8, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

colors = brewermap(10,'Set1');

figure
violinplot(trueError,{'acausal FIR', 'Blackwood', 'Zrenner', 'SSPE'}, colors)
set(gca,'Fontsize',15)
axes.LineWidth = 1.25;
ylabel('Circular Deviation (deg)')