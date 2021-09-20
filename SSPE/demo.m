% this is a demo of the causalPhaseEM_MK
% first generate some data
% creating a pure sine with some amplitude noise and added pink noise at
% sampling frequency of 1000 Hz for 10 seconds

Fs = 1000;
time = 10; 
Vlo = (10).*cos(2*pi*(5).*[1/Fs:1/Fs:time]); % random movement in amplitude, creates a little bump in PSD

[pn] = make_pink_noise(1,1e4,1/Fs);
pn = 10*pn;
data = Vlo + pn;
truePhase =  wrapToPi((2*pi*(6).*[1/Fs:1/Fs:time]))';

%%

% following helps to initialize more reasonable estimates for the scaling
% parameter 'a' and the variances for the state vector
% [init_f, init_a,init_sigma,R0] = initializeParams(data,1,1,2000); % s

% setting up initial parameters to start the causalPhaseEM code
initParams.freqs = 4; % initialization using above not great at identifying starting freq
initParams.Fs = 1000;
initParams.ampVec = .99; % in a pinch this can be initialized to 0.99 to start
initParams.sigmaFreqs = 10; % its important to use a value of this that is in the same ballpark scale
initParams.sigmaObs = 1;
initParams.window = 2000;
initParams.lowFreqBand = [4,8];

[phase,phaseBounds, fullX] = causalPhaseEM_MKmdl(data, initParams);
phase = reshape(phase', size(phase,1) * size(phase,2),1);
phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));

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
title(['Causal Phase Estimate with SP-EM, Error = ',num2str(err_Spcausal)])
xlim([2000,10000])
