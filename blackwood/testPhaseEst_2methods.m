% sim2
% Here we show that when there are two oscillations in the low freqeucny
% band we struggle to estimate the right phase, however the SP/MK approach
% does great at pulling this out.


% lets test both the hilbert transformer and the current method jointly:
% forcing data sent to both to be the same.

Fs  = 1000;
time= 15;
freq = 6;
cnt = 1;
mult = [0:4,5,4:-1:0];

ang_var2dev = @(v)(sqrt(-2*log(1-v)));

timeDelayAWMeth = zeros (10,11);
% timeDelayMKMeth = zeros(10,11);
circVarAWMeth = zeros(10,11);
% circVarMKMeth = zeros(10,11);

for secOscDist = [-5:5]
    tic
    ampOfSecOsc = [0.25:.5:5];
    parfor cnt2 = 1:length(ampOfSecOsc)      
        initParams = struct();
        data = (10 ).*cos(2*pi*(freq)*[1/Fs:1/Fs:time]) +...
            ampOfSecOsc(cnt2)*(10).*cos(2*pi*(freq+secOscDist)*[1/Fs:1/Fs:time] + pi/3) ...
            + .25*randn(1,time*Fs);
        realPhase =wrapToPi(2*pi*(freq)*[1/Fs:1/Fs:time]);
        
%         [init_f, init_a,init_sigma,R0] = initializeParams(data,2,1,5000);

        initParams.freqs = [4+(.17*(mult(cnt))),8-(.17*(mult(cnt)))];
        initParams.Fs = 1000;
        initParams.ampVec = [.99,.99];
        initParams.sigmaFreqs = init_sigma;
        initParams.sigmaObs = 1;
        initParams.window = 5000;
        initParams.lowFreqBand = [4+(.17*(mult(cnt))),8-(.17*(mult(cnt)))];

%         [phase,phaseBounds, fullX] = causalPhaseEM_MKmdl(data, initParams);

%         phase = reshape(phase', size(phase,1) * size(phase,2),1);
%         phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));
%         fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));
% 
        [phase_new,lp,up, analytic] = hilbert_transformer_phase(data', 183);
                   
        timeDelayAWMeth(cnt2,cnt) = angle(mean(exp(1i*(phase_new(5001:end)' - realPhase(5001:end)))))*(180/pi);
%         timeDelayMKMeth(cnt2,cnt) = angle(mean(exp(1i*(phase(5001:end)' - realPhase(5001:end)))))*(180/pi);

        % circular variance
        circVarAWMeth(cnt2,cnt) =1-abs(mean(exp(1i*(phase_new(5001:end)' - realPhase(5001:end)))));
%         circVarMKMeth(cnt2,cnt) =1-abs(mean(exp(1i*(phase(5001:end)' - realPhase(5001:end)))));
    end
    toc
    cnt = cnt + 1;
end

circVarAWMeth = rad2deg(ang_var2dev(circVarAWMeth));
% circVarMKMeth = rad2deg(ang_var2dev(circVarMKMeth));

%%
subplot(221)
imagesc(10*[0.25:.5:5], 6+[-5:5],timeDelayMKMeth(1:10,:)')
c = colorbar;
c.Label.String = 'Average Phase Error';
set(gca,'Fontsize', 14)
xlabel('Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
title('Causal Phase Estimate with MK-EM','Fontsize', 16)
caxis([0,45])

subplot(222)
imagesc(10*[0.25:.5:5], 6+[-5:5],timeDelayAWMeth(1:10,:)')
caxis([0,45])
set(gca,'Fontsize', 14)
xlabel('Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
c = colorbar;
c.Label.String = 'Average Phase Error';
title('Anderson-Widge method', 'Fontsize', 16)

subplot(223)
imagesc(10*[0.25:.5:5], 6+[-5:5],circVarMKMeth(1:10,:)')
caxis([0,180])
set(gca,'Fontsize', 14)
xlabel('Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
c = colorbar;
c.Label.String = 'Average Phase Error';
 title('Circular Variance', 'Fontsize', 16)

subplot(224)

imagesc(10*[0.25:.5:5], 6+[-5:5],circVarAWMeth(1:10,:)')
caxis([0,180])
set(gca,'Fontsize', 14)
xlabel('Amplitude', 'Fontsize', 18)
ylabel('Frequency (Hz)', 'Fontsize', 18)
c = colorbar;
c.Label.String = 'Average Phase Error';
 title('Error in Phase (degrees)', 'Fontsize', 16)

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
