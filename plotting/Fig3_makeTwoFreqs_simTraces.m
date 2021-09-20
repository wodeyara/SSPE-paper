secOscDist = -1;
ampOfSecOsc = [1.8];
freq = 6; 
Fs = 1000;
time = 15;

initParams = struct();
data = (25).*cos(2*pi*(freq)*[1/Fs:1/Fs:time]) +...
ampOfSecOsc*(25).*cos(2*pi*(freq+secOscDist)*[1/Fs:1/Fs:time] + pi/3) ...
+ .5*randn(1,time*Fs);

truePhase =wrapToPi(2*pi*(freq)*[1/Fs:1/Fs:time]);
truePhaseSecOsc  = wrapToPi(2*pi*(freq+secOscDist)*[1/Fs:1/Fs:time] + pi/3);

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

figure
% subplot(311)
% subplot(312)
plot(truePhase(7001:10000), 'Linewidth',3)
hold on
plot(truePhaseSecOsc(7001:10000), 'Linewidth',3)
figure
% subplot(313)
plot(truePhase(7001:8000), 'Linewidth',3)
hold on
plot(filt_phase(7001:8000), 'Linewidth',1.5)
plot(phase_AW(7001:8000), 'Linewidth',1.5)
plot(zrenner_phase([7001:8000]-999), 'Linewidth',1.5)
plot(phase(7001:8000), 'Linewidth',2, 'Linestyle', '--', 'Color', 'red')
legend({'True','acausal FIR','Blackwood','Zrenner','SSPE'})
set(gca,'Fontsize', 16)
ylabel('Phase')
xlabel('Time (ms)')

%%
allColors = brewermap(5,'Set1');
figure
plot(truePhase(7001:10000), 'Linewidth',3,'color',allColors(1,:))
hold on
plot(truePhaseSecOsc(7001:10000), 'Linewidth',3,'color',allColors(2,:))
ylabel('Phase')

yyaxis right
plot(data(7000:10000),'Linewidth',3,'color','black')
ylabel('Signal')
xlabel('Time')