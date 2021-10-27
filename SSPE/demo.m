% this is a demo of the causalPhaseEM_MK
% first generate some data
% creating a pure sine and added pink noise at
% sampling frequency of 1000 Hz for 10 seconds

Fs = 1000;
time = 10; 
Vlo = (10).*cos(2*pi*(5).*[1/Fs:1/Fs:time]); % random movement in amplitude, creates a little bump in PSD

[pn] = make_pink_noise(1,1e4,1/Fs);
pn = 10*pn;
data = Vlo + pn;
truePhase =  wrapToPi((2*pi*(6).*[1/Fs:1/Fs:time]))';

%%

% setting up initial parameters to start the causalPhaseEM code
initParams.freqs = 4; % only tracking a single oscillator
initParams.Fs = 1000; % sampling frequency for data
initParams.ampVec = .99; % in a pinch this can be initialized to 0.99 to start - represents rate of damping
initParams.sigmaFreqs = 10; % its important to use a value of this that is in the same ballpark scale, estimate of covariance
initParams.sigmaObs = 1; % observation variance for measurement noise
initParams.window = 2000; % the analysis uses 2000 samples from start of data to fit params
initParams.lowFreqBand = [4,8]; % band that we expect oscillator of interest to be in

[phase,phaseBounds, fullX,returnParams] = causalPhaseEM_MKmdl_noSeg(data, initParams,0);
disp('Estimates of the parameters are below:')
returnParams
%%
figure
err_Spcausal = angle(mean(exp(1i*(phase(2001:end) - truePhase(2001:end)))))*(180/pi);
imagesc( 1:10000,-3:3, (abs(fullX(:,1) + 1i*fullX(:,2)))', 'AlphaData', .5)
colormap(summer)
hold on; plot(1:10000, phase, 'Linewidth', 2, 'Color', 'b')
plot(squeeze(phaseBounds(:,1)), 'Linewidth', 2, 'color','r')
hold on
plot(squeeze(phaseBounds(:,2)), 'Linewidth', 2,'color','r')
set(gca,'Fontsize', 16)
xlabel('Time (in samples)')
ylabel('Phase')
h = colorbar;
ylabel(h,'Amplitude')
title(['Causal Phase Estimate with SSPE, Error = ',num2str(err_Spcausal)])
xlim([2000,10000])
