% This is example code that analyses a vector of data. It was set up for
% LFP data with 250 s of theta and the code sets out to track that. 

%%
% re-estimate parameters for every 10 seconds
% split data into 40 segments of 10 seconds separated by 1 second every
% time and then do a parfor across all of them, then do next 40 and so on.
% Save all initParams and apply them in a separate function to data. 
%% setting up data and getting SNR
%basic parameters:
Fs = 1000;
time = 10;
thetaBand = [8,14];

% creating the matrix that will be fit 
data = allData{5};
data = data(1:300e3);
dataToFit = zeros(length(data)/Fs - 9,Fs*time);
for i = 1:length(data)/Fs - 9
    dataToFit(i,:) = data((i-1)*Fs + 1 : (i-1)*Fs + Fs*time);
end

% calculate SNR for each 10 second Interval
peak_SNR = zeros(length(data)/Fs - 9,1);
for i = 1:length(data)/Fs - 9
    epochs = reshape(dataToFit(i,:)', Fs, time);
    [peak_SNR(i), peak_frequency] = estimate_SNR_theta(epochs, Fs, ...
                                thetaBand, Fs, 0);
end

%% fit the model on all 300 segments and save out the initial parameters. 

allInitParams = {};
parfor i = 1:length(data)/Fs - 9
tic
    initParams = struct();
    initParams.Fs = 1000;
    initParams.lowFreqBand = thetaBand;
    initParams.window = Fs*time;
    [initParams.freqs, initParams.ampVec, initParams.sigmaFreqs, initParams.sigmaObs, ...
        stateVec, stateCov] = fit_MKModel_multSines(squeeze(dataToFit(i,:)),[1, 7, 60], ...
                                                    1000,[.99, 0.99,.99], [20, 20,20,5],1);
    allInitParams{i} = initParams;
    allStateVec{i} = stateVec;
toc
end
save('appToData_2.mat', 'allInitParams','allStateVec','dataToFit','data','peak_SNR')
%% apply the model
% only use parameters estimated before the 10 second before the second to
% which we are applying parameters so for 51st second only use 40 - 50  and everthing before so
% we never have the second itself incorporated into parameter estimation.
secsToUse = 50:(length(data)/Fs) -1;
for i = 1:length(secsToUse)
    tic
    for j = 1:secsToUse(i)-9
        % for 50 we go up till index 41 which has 40001:50000
        % and we will be using 50001:51000 as the epoch and so on...
        initParams = allInitParams{j};
        data_onesec = data(secsToUse(i)*Fs + 1: secsToUse(i)*Fs + Fs);
        [~, netError(i,j)] = runSSPE_KalmanFilter(data_onesec,initParams,0); 
        sumSSE = sum(data_onesec.^2);
        normNetError(i,j) = netError(i,j)/sumSSE;
    end
    toc
end


%% identify moments when we found theta
wasthereLowFreq = zeros(length(allInitParams),1);
for i = 1:length(allInitParams)
    tmp = allInitParams{i};
   wasthereLowFreq(i) = any(tmp.freqs> tmp.lowFreqBand(1) & tmp.freqs< tmp.lowFreqBand(2));
end
save('appToData_2.mat', 'netError', 'normNetError', 'wasthereLowFreq', '-append')