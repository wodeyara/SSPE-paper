function [phase,phaseBounds,allX_full,phaseWidth,returnParams] = causalPhaseEM_MKmdl(y,initParams)
% EXISTS FOR BACKWARD COMPATIBILITY -> USE causalPhaseEM_MKmdl_noSeg()
% causal phase estimates using the SP model and EM with fixed interval
% smoothing across windows of data
% primary benefit: assume a transitory burst of oscillatory activity in the
% range of your bandpass filter/ assume a peak shift in the data towards the edge
% of the band pass filter. These are problems unaddressed for instantaneous
% phase estimation right now

% Algorithm:
% after estimating reasonable initialization points we need to run EM on
% the data - at whatever rate makes it possible to run it again before
% getting the next window of data. So while it might be slow here, a C++
% implementation is likely going to be much faster meaning we have have
% small windows (up to the frequency limits of course)
% INPUT: 
% y - data/observation
% initParams - a structure containing the following parameters: 
%           .freqs - oscillator center frequencies to be tracked
%           .Fs - sampling frequency
%           .ampVec - decay rate for each oscillator (initialize to 0.99)
%           .sigmaFreqs - variance for each oscillator (.1 is good default)
%           .sigmaObs - variance for observation (1)
%           .windowSize - how large is the window used to fit the parameters
%           .lowFreqBand - where is the oscillator you wish to track, give
%           a frequency range
% OUTPUT
% phase - Estimated phase for the first oscillator in the frequency band range 
% phaseBounds - The phase credible intervals, use to judge confidence in the phase
% allX_full - Estimated state for the first oscillator in freq band range
% phaseWidth - Single value in degrees informing us of confidence in phase
% returnParams - Estimated parameter values in same format at initParams
% Last edit: Ani Wodeyar 6/29/2021
freqs = initParams.freqs;
Fs = initParams.Fs;
ampVec = initParams.ampVec;
sigmaFreqs = initParams.sigmaFreqs;
sigmaObs = initParams.sigmaObs;
windowSize = initParams.window;
lowFreqBand = initParams.lowFreqBand;

if windowSize < Fs
    disp('The window size needs to be different. Setting it equal to sampling rate')
    windowSize = Fs;
end

if length(y) < 2*windowSize
    disp('please enter a larger vector of observations (should be at leas 2x as big as window size)')
    return
end

numSegments = floor(length(y)/windowSize);
ang_var2dev = @(v) sqrt(-2*log(v)); % note the difference in definition (ie not (1-v))

data = y(1:windowSize);
% first run to set up parameters
[omega, ampEst, allQ, R, stateVec, stateCov] = fit_MKModel_multSines(data,freqs, Fs,ampVec, sigmaFreqs,sigmaObs);
lowFreqLoc = find((omega>lowFreqBand(1)) & (omega<lowFreqBand(2)),1);
returnParams.freqs = omega;
returnParams.ampVec = ampEst;
returnParams.sigmaFreqs = allQ;
returnParams.sigmaObs = R;

 
if isempty(lowFreqLoc)
    disp('Low freq band limits incorrect OR there is no low freq signal; retaining initial params')
    omega = freqs;
    ampEst = ampVec;
    allQ = sigmaFreqs;
    [~,lowFreqLoc] = min(abs(freqs-mean(lowFreqBand))); % pick frequency closest to middle of low frequency range
end

% for loop that runs through rest of the data reestimating parameters after
% generating phase estimates for the whole period using past parameter ests
% and the kalman filter
[phi, Q, M] = genParametersSoulatMdl_sspp(omega, Fs, ampEst, allQ);
phase = zeros(numSegments, windowSize);
phaseBounds = zeros(numSegments, windowSize,2);
allX_full = zeros(numSegments, windowSize, 2);
phaseWidth = zeros(numSegments,windowSize);

for seg = 2:numSegments
%     tic
    y_thisRun = y((seg-1)*windowSize + 1: seg*windowSize);
    % running Kalman filter over one window before re-running EM
    % start below with the end of the EM run x and stae cov
    allX = zeros(length(freqs)*2, windowSize);
    allP = zeros(length(freqs)*2,length(freqs)*2, windowSize);
    
    x = stateVec(:,end);
    P = squeeze(stateCov(:,:,end));
    
    for i = 1:(length(y_thisRun))
        % kalman update
        [x_new,P_new] = oneStepKFupdate_sspp(x,y_thisRun(i),phi,M,Q,R,P);
        allX(:,i) = x_new;
        P_new = (P_new + P_new') /2; % forcing symmetry to kill off rounding errors
        allP(:,:,i) = P_new; 
        
        % estimate phase
        phase(seg, i) = angle(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
        samples = mvnrnd(x_new(lowFreqLoc*2-1:lowFreqLoc*2),...
            P_new(lowFreqLoc*2-1:lowFreqLoc*2,lowFreqLoc*2-1:lowFreqLoc*2),2000);
        
        sampleAngles = (angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(seg,i)))); % removing mean
        lowerBnd = (prctile(sampleAngles,2.5));
        upperBnd = (prctile(sampleAngles,97.5));
        phaseBounds(seg,i,:) = sort([lowerBnd + (phase(seg,i)), ...
                                     upperBnd + (phase(seg,i))]); % can have a range of [0,2pi]
        phaseWidth(seg,i) = rad2deg(ang_var2dev(abs(mean(exp(1i*sampleAngles)))));

        % update state and state cov
        P = P_new;
        x = x_new;
    end
    
    allX_full(seg,:,:) = allX(lowFreqLoc*2-1:lowFreqLoc*2,:)';
    
end
    

