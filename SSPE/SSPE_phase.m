function [phase,phaseBounds,newAllX] = SSPE_phase(data,initParams,lowFreqLoc)

freqs = initParams.freqs;
Fs = initParams.Fs;
ampVec = initParams.ampVec;
sigmaFreqs = initParams.sigmaFreqs;
sigmaObs = initParams.sigmaObs;

y = data;
ang_var2dev = @(v)(sqrt(-2*log(v)));

% need to initialize all my parameters

[phi, Q, M] = genParametersSoulatMdl_sspp(freqs, Fs, ampVec, sigmaFreqs); 
R = sigmaObs;

xstart = zeros(2*length(freqs),1);
Pstart = .001 * eye(2*length(freqs));
x = xstart;
P = Pstart;

%% run forward KF
for i = 1:(length(y))
    
    [x_new,P_new] = oneStepKFupdate_sspp(x,y(i),phi,M,Q,R,P);
    allX(:,i) = x_new;
    allP(:,:,i) = P_new;
    
    P = P_new;
    x = x_new;
end

% now we have P_new as P_N_N and P as P_(N-1)_(N-1)
% same for x
% can now run the fixed interval smoother

%% run the fixed interval smoother across all data

x_n = x_new;
P_n = P_new;

newAllX(:, length(y)) = allX(:,length(y));
newAllP(:,:,length(y)) = allP(:,:,length(y));
for i = length(y)-1:-1:1
    
    [x_backone,P_backone, J_one] =  fixedIntervalSmoother_sspp(x, x_n, P, P_n, phi, Q);
    x_n = x_backone;
    P_n = triu(P_backone,1) + triu(P_backone,0)';
    x = allX(:,i);
    P = squeeze(allP(:,:,i));    
    
    newAllX(:,i+1) = x_backone;
    newAllP(:,:,i+1) = P_backone;
    allJ(:,:,i+1) = J_one;
end
[x_backone,P_backone, J_one] =  fixedIntervalSmoother_sspp(x, x_n, P, P_n, phi, Q);
newAllX(:,1) = x_backone;
newAllP(:,:,1) = P_backone;
allJ(:,:,1) = J_one;

phase = zeros(length(y),length(freqs));
% pre-sample the multivariate gaussian
allRandVals = randn(length(freqs)*2, 2000);
parfor i = 1:length(y)
    x_n = newAllX(:,i);
    P_n = newAllP(:,:,i);
    
    P_n_decomp = chol(P_n);
    
    phase(i,:) = angle(x_n(1:2:end) + 1i*x_n(2:2:end));
    samples = x_n + P_n_decomp * allRandVals; % samps = mu + Cov_decomp*(unit normal samps);
%     samples = mvnrnd(x_n,P_n,2000); 
    tmpPhase = angle(x_n(lowFreqLoc*2-1) + 1i*x_n(lowFreqLoc*2));
    sampleAngles = (angle(exp(1i*angle(samples(lowFreqLoc*2-1,:) + ...
                            1i*samples(lowFreqLoc*2,:)) - 1i*tmpPhase))); % removing mean
    
%     lowerBnd = (prctile(sampleAngles,2.5));
    sampMean = mean(sampleAngles);
    sampSTD = 2*std(sampleAngles);
    lowerBnd = sampMean - sampSTD;
    upperBnd = sampMean + sampSTD;
    
    phaseBounds(i,:) = sort([lowerBnd + (tmpPhase), ...
                                 upperBnd + (tmpPhase)]); % can have a range of [0,2pi]

    
end


% need to estimate P_t_(t-1) for t = 1:N
P_tmp = phi * squeeze(allP(:,:,end-1)) * phi' +Q;
K = P_tmp * M * (1/(M' * P_tmp*M + R));


