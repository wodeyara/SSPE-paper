% creating a fake signal to test the hilbert_Transformer_phase code
clear all 
close all

Fs  = 3e4;
time= 10;
freq = 6;
cnt = 1;
pinkNoise = dsp.ColoredNoise('pink',Fs,10);
timeDelay = [];

for secOscDist = -5.5:.5:5.5
    cnt2 = 1;
    tic
    for ampOfSecOsc = 0.5:.5:10       
        data = cos(2*pi*(freq)*[1/Fs:1/Fs:time]) + ampOfSecOsc*cos(2*pi*(freq+secOscDist)*[1/Fs:1/Fs:time] + pi/4) ...
            + reshape(step(pinkNoise),1,Fs*time);
        realPhase = wrapToPi((2*pi*(freq)*[1/Fs:1/Fs:time]));
        
        
        
        [phase, estimate_mask, analytic] = hilbert_transformer_phase(data, 10000);
        % getting metrics of interest
        dsRealPhase = wrapToPi((2*pi*(freq)*[1/500:1/500:time])); %500 is downsampled sampling frequency
        timeDelay(cnt2,cnt) = (mean((abs(phase - dsRealPhase'))/(2*pi) * (1/freq) *1000));

        cnt2 = cnt2 + 1;
    end
    toc
    cnt = cnt + 1;
end
%%
dsRealPhase = wrapToPi((2*pi*(freq)*[1/500:1/500:time]));
subplot(121)
plot(dsRealPhase)
hold on
plot(phase)
legend({'Orig Phase', 'Est Phase'})
subplot(122)
rose((phase - dsRealPhase'),30)
title(num2str(mean((abs(phase - dsRealPhase'))/(2*pi) * (1/freq) *1000)) ,'Fontsize', 18)

%% what happens when the signals aren't sinusoidal? (but still biologically feasible - how to simulate?)


%% COMMENTS
% under white noise, the system works great even under considerable noise
% (STD of 20 for white noise was the limit I checked) . Performance falls
% apart of course when the frequency of interest is 3 Hz. 

% if we have a second oscillation which is at the edge of the pass band -
% 9/10 Hz then the performance falls spectacularly. At 10 Hz however,
% performance falls when the amplitude of that oscillation is several times
% higher than the oscillation of interest.
