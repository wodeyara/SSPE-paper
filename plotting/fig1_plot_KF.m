

function [h1,h2] = plot_KF(y,initParams)

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
%%
plot(1/Fs:1/Fs:3,y(1:3000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
hold on
plot(1/Fs:1/Fs:length(data)/Fs,stateVec(1,:),'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',2)
set(gca,'YTickLabel', {},'XTickLabel', {})

%%
% for loop that runs through rest of the data reestimating parameters after
% generating phase estimates for the whole period using past parameter ests
% and the kalman filter
[phi, Q, M] = genParametersSoulatMdl_sspp(omega, Fs, ampEst, allQ);
phase = zeros(numSegments, windowSize);
phaseBounds = zeros(numSegments, windowSize,2);
allX_full = zeros(numSegments, windowSize, 2);
circstd = zeros(numSegments,windowSize);
cosineSum = zeros(numSegments,windowSize);

for seg = 2:numSegments
%     tic
    y_thisRun = y((seg-1)*windowSize + 1: seg*windowSize);
    % running Kalman filter over one window before re-running EM
    % start below with the end of the EM run x and stae cov
    allX = zeros(length(freqs)*2, windowSize);
    allP = zeros(length(freqs)*2,length(freqs)*2, windowSize);
    
    x = stateVec(:,end);
    P = squeeze(stateCov(:,:,end));
    
    for i = 1:(length(y_thisRun)-1)
        % kalman update
        [x_new,P_new] = oneStepKFupdate_sspp(x,y_thisRun(i),phi,M,Q,R,P);
        allX(:,i) = x_new;
        P_new = (P_new + P_new') /2; % forcing symmetry to kill off rounding errors
        allP(:,:,i) = P_new; 
        
        
        % estimate phase
        phase(seg, i) = angle(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
        samples = mvnrnd(x_new(lowFreqLoc*2-1:lowFreqLoc*2),...
            P_new(lowFreqLoc*2-1:lowFreqLoc*2,lowFreqLoc*2-1:lowFreqLoc*2),2000);
        % akin to what is done under Technique 2 on page 206 of FIsher 1993
        % statistical analysis of circular data
        sampleAngles = angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(seg,i))); % removing mean
        tmpAngleStd = abs(wrapToPi((prctile(sampleAngles,97.5) - prctile(sampleAngles,2.5)))/2); %[0,pi/2]
        phaseBounds(seg,i,:) = sort([(phase(seg,i)- tmpAngleStd), (phase(seg,i) + tmpAngleStd)]); % can have a range of [0,pi]
        circstd(seg,i) = ang_var2dev(abs(mean(exp(1i*sampleAngles))));
        cosineSum(seg,i) = mean(cos(2*(sampleAngles)));
%         sampleAmp = abs(samples(:,1) + 1i*samples(:,2));
%         amp(seg,i) = std(sampleAmp);

        % update state and state cov
        P = P_new;
        x = x_new;
    end
    % ending it on N-1
    [x_new,P_new] = oneStepKFupdate_sspp(x,y_thisRun(i),phi,M,Q,R,P);
    allX(:,i+1) = x_new;
    allP(:,:,i+1) = P_new;
    
    % estimate phase
    phase(seg, i) = angle(x_new(lowFreqLoc*2-1) + 1i* x_new(lowFreqLoc*2));
    samples = mvnrnd(x_new(lowFreqLoc*2-1:lowFreqLoc*2), P_new(lowFreqLoc*2-1:lowFreqLoc*2,lowFreqLoc*2-1:lowFreqLoc*2),2000);
    sampleAngles = angle(exp(1i*angle(samples(:,1) + 1i*samples(:,2)) - 1i*phase(seg,i)));
    tmpAngleStd = abs(wrapToPi((prctile(sampleAngles,97.5) - prctile(sampleAngles,2.5)))/2);
    phaseBounds(seg,i,:) = sort([(phase(seg,i)+ tmpAngleStd), (phase(seg,i) + tmpAngleStd)]);
    circstd(seg,i) = ang_var2dev(abs(mean(exp(1i*sampleAngles))));%exp(1i*sampleAngles)
    cosineSum(seg,i) = mean(cos(2*(sampleAngles)));
    

    allX_full(seg,:,:) = allX(lowFreqLoc*2-1:lowFreqLoc*2,:)';
    %%
    plot(2+(1/Fs):1/Fs:3,y_thisRun(1:1000), 'Color',[0.8500, 0.3250, 0.0980], 'Linewidth',2);
    hold on
    plot(2+(1/Fs):1/Fs:3,allX(1,1:1000), 'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',3,'LineStyle','-.')
    plot(2+(1/Fs):1/Fs:3,allX(2,1:1000), 'Color',[34,139,34]/255, 'Linewidth',3,'LineStyle','-.')
    grid on
    set(gca,'YTickLabel', {},'XTickLabel', {})
    legend({'Observation', 'Real','Imaginary'})
    set(gca,'Fontsize', 16)
    hAxis = gca
    hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
    hAxis.XRuler.FirstCrossoverValue = 0; % X crossover with Z axis
    hAxis.YRuler.SecondCrossoverValue = 0; % X crossover with Z axis
    hAxis.YRuler.FirstCrossoverValue = 0; % X crossover with Z axis
    hAxis.XTick = 2:.1:3
    hAxis.TickDir = 'both';
    hAxis.GridColor = 	[0, 0.4470, 0.7410];
%%

plot(2+(1/Fs):1/Fs:3, phase(2,1:1000),'Linewidth',2);
ylabel('Phase')
ylim([-pi,pi])
h = gca
set(h,'Fontsize', 16);
hold on
plot(2+(1/Fs):1/Fs:3, squeeze(phaseBounds(2,1:1000,1)),'Linewidth',2,'Linestyle','--','Color','red')
plot(2+(1/Fs):1/Fs:3, squeeze(phaseBounds(2,1:1000,2)),'Linewidth',2,'Linestyle','--','Color','red')
ylim([-pi*1.2,pi*1.2])
h.XAxis.Visible = 'off'
%%
    % update below to take in the updated x start and cov start 
   [freqs, ampVec, sigmaFreqs, R, stateVec, stateCov,] = fit_MKModel_multSines(y_thisRun,omega, Fs, ampEst, allQ, R);
   tmp = find((freqs>lowFreqBand(1)) & (freqs<lowFreqBand(2)),1);
   
    if isempty(tmp)
        disp('Low freq band limits incorrect OR there is no low freq signal; retaining old parameters')       
    else
        lowFreqLoc = tmp;
        omega = freqs;
        ampEst = ampVec;
        allQ = sigmaFreqs;
    end
   
   [phi, Q, M] = genParametersSoulatMdl_sspp(omega, Fs, ampEst, allQ);
%     toc
end
    

