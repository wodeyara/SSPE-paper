function [phase, netError,phaseBounds,phaseWidth, allX,amp] = runSSPE_KalmanFilter(data,initParams,credIntflag)
% given the parameters in initParams, this will run the kalman filter over
% data using the given parameters
freqs = initParams.freqs;
Fs = initParams.Fs;
ampVec = initParams.ampVec;
sigmaFreqs = initParams.sigmaFreqs;
sigmaObs = initParams.sigmaObs;

ang_var2dev = @(v)(sqrt(-2*log(v)));

lowFreqLoc = find(freqs>initParams.lowFreqBand(1) & freqs<initParams.lowFreqBand(2));
if isempty(lowFreqLoc)
    disp('setting arbitrary low freq location, dont trust phase')
    lowFreqLoc = 1;
end

[phi, Q, M] = genParametersSoulatMdl_sspp(freqs, Fs, ampVec, sigmaFreqs);

x = zeros(length(freqs)*2,1);
P = eye(length(freqs)*2);

allX = zeros(length(freqs)*2, Fs);
allP = zeros(length(freqs)*2,length(freqs)*2, Fs);
phaseBounds = [];
phaseWidth = [];
for i = 1:(length(data))
    % kalman update
    [x_new, P_new] = oneStepKFupdate_sspp(x,data(i),phi,M,Q,sigmaObs,P);
    allX(:,i) = x_new;
    P_new = (P_new + P_new') /2; % forcing symmetry to kill off rounding errors
    allP(:,:,i) = P_new; 
    
    % estimate phase
    amp(i) = abs(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
    phase(i) = angle(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
    if credIntflag ==1
    samples = mvnrnd(x_new(lowFreqLoc*2-1:lowFreqLoc*2),...
        P_new(lowFreqLoc*2-1:lowFreqLoc*2,lowFreqLoc*2-1:lowFreqLoc*2),2000);
    sampleAngles = (angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(i)))); % removing mean
    lowerBnd = (prctile(sampleAngles,2.5));
    upperBnd = (prctile(sampleAngles,97.5));
    phaseBounds(i,:) = sort([lowerBnd + (phase(i)), ...
                                 upperBnd + (phase(i))]); % can have a range of [0,2pi]
    phaseWidth(i) = rad2deg(ang_var2dev(abs(mean(exp(1i*sampleAngles)))));
    end

    % update state and state cov
    P = P_new;
    x = x_new;
end
% ending it on N-1
% [x_new,P_new] = oneStepKFupdate_sspp(x,data(i),phi,M,Q,sigmaObs,P);
% allX(:,i+1) = x_new;
% allP(:,:,i+1) = P_new;

% estimate phase
% phase(seg, i) = angle(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
% samples = mvnrnd(x_new(lowFreqLoc*2-1:lowFreqLoc*2), P_new(lowFreqLoc*2-1:lowFreqLoc*2,lowFreqLoc*2-1:lowFreqLoc*2),2000);
% sampleAngles = angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(seg,i)));
% tmpAngleStd = abs(wrapToPi((prctile(sampleAngles,97.5) - prctile(sampleAngles,2.5)))/2);
% phaseBounds(seg,i,:) = sort([(phase(seg,i)+ tmpAngleStd), (phase(seg,i) + tmpAngleStd)]);
netError = [];
% pred_Y = sum(allX(1:2:end,:),1);
% netError = sum((data-pred_Y).^2);
