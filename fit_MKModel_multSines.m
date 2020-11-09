
function [omega, ampEst, allQ, R,stateVec, stateCov] = fit_MKModel_multSines(data,freqs, Fs,ampVec, sigmaFreqs,sigmaObs)
% fitting the soulat model using the shumway-stoffer method (basically
% assuming no harmonics); EM has an analytic solution that needs the
% covariance structure between adjacent time points (this is what we need
% the fixed interval smoother for) 
% % This is the extension of the fitSoulatModel code to fit multiple
% sinusoids, still with no harmonics. 
% INPUTS:
% data - currently only accepts vector data so n x 1
% freqs - initialize with some frequencies that are appropriate for data
% Fs - sampling frequency
% sigmaFreqs - what are variances at each of the frequencies specified
% sibmaObs - what is the expected observation noise variance?
% Ani Wodeyar 3/17/2020
% have now fixed the indexing of the fixed interval smoother, was off by 1
% before
% Ani Wodeyar 9/23/2020

if isempty(sigmaFreqs)
    sigmaFreqs = 0.1*ones(length(freqs),1);
end
if isempty(sigmaObs)
    sigmaObs = 10;
end
y = data;

freqEst = freqs/Fs;
% need to initialize all my parameters

[phi, Q, M] = genParametersSoulatMdl_sspp(freqs, Fs, ampVec, sigmaFreqs); 
R = sigmaObs;

iter = 1;
errorVal = Inf;
%iterate though maximizing likelihood

while iter < 400 && errorVal(iter)>1e-3

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
    P_n = P_backone;
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

% need to estimate P_t_(t-1) for t = 1:N
P_tmp = phi * squeeze(allP(:,:,end-1)) * phi' +Q;
K = P_tmp * M * (1/(M' * P_tmp*M + R));
P_N_N1 = (eye(length(P_tmp)) - K * M') * phi * squeeze(allP(:,:,end-1));

allP_N_N1(:,:,length(y)) = P_N_N1;

for i = length(y)-1:-1:2
   allP_N_N1(:,:,i)=squeeze(allP(:,:,i)) * squeeze(allJ(:,:,i-1))' + ...
                    squeeze(allJ(:,:,i))*(squeeze(allP_N_N1(:,:,i+1)) - phi * squeeze(allP(:,:,i))) * squeeze(allJ(:,:,i-1))';
end
allP_N_N1(:,:,1) = eye(size(squeeze(allP(:,:,1))));

%% update the parameter estimate from optimizing exp cond likelihood
% the portion here is a direct lift from shumway and stoffer rather than
% Soulat and Purdon

A = Pstart + xstart * xstart';
B = zeros(size(newAllP,1),size(newAllP,2));
C = zeros(size(newAllP,1),size(newAllP,2));
R = zeros(1);

for i = 1:length(y)
    if i > 1
        A = A + squeeze(newAllP(:,:,i-1)) + newAllX(:,i-1) *newAllX(:,i-1)';
        B = B + squeeze(allP_N_N1(:,:,i)) + newAllX(:,i) *newAllX(:,i-1)';
    end
    C = C + squeeze(newAllP(:,:,i)) + newAllX(:,i) *newAllX(:,i)';
    R = R+ M'* squeeze(newAllP(:,:,i)) * M + (y(i) - M'*newAllX(:,i)) *(y(i) - M'*newAllX(:,i))' ;
end

R = (1/(length(y))) * (R);

oldFreq = freqEst * 1000 /(2*pi);

freqEst = zeros(1,length(freqs));
ampEst = zeros(1,length(freqs));
allQ = zeros(1,length(freqs));

for numFreqs = 1: length(freqs)
    B_tmp = B((numFreqs-1)*2 + 1: numFreqs*2,(numFreqs-1)*2 + 1: numFreqs*2);
    A_tmp = A((numFreqs-1)*2 + 1: numFreqs*2,(numFreqs-1)*2 + 1: numFreqs*2);
    C_tmp = C((numFreqs-1)*2 + 1: numFreqs*2,(numFreqs-1)*2 + 1: numFreqs*2);
    freqEst(numFreqs) = (atan((B_tmp(2,1) - B_tmp(1,2))/trace(B_tmp)));
    ampEst(numFreqs) = min(sqrt((B_tmp(2,1) - B_tmp(1,2))^2 + trace(B_tmp)^2)/trace(A_tmp),1-eps); % constrained to 1
    allQ(numFreqs) = 1/(2*length(y)) * (trace(C_tmp) - ampEst(numFreqs).^2 * trace(A_tmp));
end
    
[phi, Q] = genParametersSoulatMdl_sspp(freqEst * 1000 /(2*pi), Fs, ampEst, allQ); 

omega = freqEst * 1000 /(2*pi);
stateVec = newAllX;
stateCov = newAllP;

% xstart = newAllX(:,1);
% Pstart = squeeze(newAllP(:,:,1));

iter = iter + 1;
errorVal(iter) =sum(abs(omega - oldFreq));
end

