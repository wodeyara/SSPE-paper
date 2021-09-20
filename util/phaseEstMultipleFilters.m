
%% define variables and functions
fs = 1000;
peak_frequency = 7;

ang_var2dev = @(v) sqrt(-2*log(1-v));
%% create filters
filter_objects = {};
        
for ord = [2 3 4 5 6 7 8] % FIR - windowed sinc
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'CutoffFrequency1', peak_frequency-1, 'CutoffFrequency2', peak_frequency+1, 'SampleRate', fs, 'DesignMethod', 'window')};
end
for ord = [3 4 5 6 7 8] % FIR - least squares (equiripple is similar)
    filter_objects = {filter_objects{:} designfilt('bandpassfir', 'FilterOrder', round(ord * (fs/peak_frequency)), 'StopbandFrequency1', peak_frequency-2, 'PassbandFrequency1', peak_frequency-1, 'PassbandFrequency2', peak_frequency+1, 'StopbandFrequency2', peak_frequency+2, 'SampleRate', fs, 'DesignMethod', 'ls')};
end
    
%% run filters

analyticSigCell = cellfun(@(x) hilbert(filtfilt(x, data)), filter_objects, 'UniformOutput', false);
for i = 1:length(analyticSigCell)
    analyticSig(i,:) = analyticSigCell{i};
end

bestPhaseEst = angle(mean(analyticSig,1));
phaseErrorEst = rad2deg(ang_var2dev(1 - abs(mean(analyticSig./abs(analyticSig),1))));
