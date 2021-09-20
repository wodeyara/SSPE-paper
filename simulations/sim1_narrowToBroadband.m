% setting up different simulated forms of data for the tests:

function testCausalPhaseEst(type)

% type = 'pink_noise';
time = 10;
Fs = 1000;
PAC = 0;
AAC = 0;

initParams.freqs = [4]; 
initParams.Fs = 1000;
initParams.ampVec = [.975];
initParams.sigmaFreqs = [1];
initParams.sigmaObs = .1;
initParams.window = 2000;
initParams.lowFreqBand = [4,8];

% now running FIR + Hilbert
fNQ = Fs/2;
locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

%determine default parameters (median of optimized values)
default_parameters = [];
default_parameters.window_length = 750;%ceil(median(T.optim_window));
default_parameters.filter_order = 192;%ceil(median(T.optim_filter_ord));
default_parameters.edge = 64;% ceil(median(T.optim_edge));
default_parameters.ar_order = 30;% ceil(median(T.optim_ar_ord));

D = designfilt('bandpassfir', 'FilterOrder', default_parameters.filter_order, ...
    'CutoffFrequency1', 4, 'CutoffFrequency2', 8, 'SampleRate', Fs, 'DesignMethod', 'window');



parfor iter = 1:1000
    tic
   [data,truePhase] = generateData(type, Fs,time, PAC, AAC);
   simData(iter).origData = data;
   simData(iter).truePhase = truePhase;
   
   % first run causalEM
    try
    [SP_phase,phaseBounds,fullX] = causalPhaseEM_MKmdl(data, initParams);
    SP_phase = reshape(SP_phase', size(SP_phase,1) * size(SP_phase,2),1);
    phaseBounds = reshape(permute(phaseBounds,[2,1,3]), size(phaseBounds,1) * size(phaseBounds,2),size(phaseBounds,3));    
    fullX = reshape(permute(fullX,[2,1,3]), size(fullX,1) * size(fullX,2),size(fullX,3));
    
    simData(iter).SP_phase = SP_phase;
    simData(iter).SP_phaseBounds = phaseBounds;
    simData(iter).SP_data = fullX;
    simData(iter).SP_failed=0;
    err_SPcausal(iter) = mean(abs(SP_phase(2001:end) - truePhase(2001:end)))/(2*pi);
    catch
        simData(iter).SP_failed=1;
    end
    
    lowAct = filtfilt(filtwts,1,data);    
    lowAct_analytic = hilbert(lowAct);
    simData(iter).filterData = lowAct_analytic;
    simData(iter).filt_phase = angle(lowAct_analytic);

    err_FIR_Hilb(iter) = mean(abs(angle(lowAct_analytic(2001:end))' - truePhase(2001:end)))/(2*pi);
    
    % now running Anderson-Widge method
    [phase_AW,~,~, analytic] = hilbert_transformer_causal(data', 1000, [4,8]);
    simData(iter).AW_phase = phase_AW;
    simData(iter).AW_data = analytic;
    
    err_AW(iter) = mean(abs(phase_AW(2001:end)' - truePhase(2001:end)))/(2*pi);
    
    epochs = create_epochs_overlapping(data,Fs);
    epochwindowmask = ((-default_parameters.window_length+1):0)+ceil(size(epochs,1)/2);

    [zrenner_phase, estamp] = phastimate(epochs(epochwindowmask,:), D, default_parameters.edge, default_parameters.ar_order, 128);
    simData(iter).zrenner_phase = zrenner_phase;
    simData(iter).zrenner_amp = estamp;
    
    err_zrenner(iter) = mean(abs(zrenner_phase(1002:end)' - truePhase(2001:9000)))/(2*pi);    
    
    var_SPcausal(iter) = 1 - abs(mean(exp(1i*(truePhase(2001:9000) - SP_phase(2001:9000)))));
    var_AW(iter) = 1 - abs(mean(exp(1i*(truePhase(2001:9000) - phase_AW(2001:9000)'))));
    var_FIR_Hilb(iter) = 1 - abs(mean(exp(1i*(truePhase(2001:9000) - angle(lowAct_analytic(2001:9000))'))));
    var_Zrenner(iter) = 1 - abs(mean(exp(1i*(truePhase(2001:9000) - zrenner_phase(1002:end)'))));

    toc       
end
   try 
    save(['causalPhaseEst_sim_',type,'.mat'],'simData', 'err_AW', 'err_FIR_Hilb', 'err_SPcausal','err_zrenner',...
                                              'var_SPcausal','var_AW','var_FIR_Hilb','var_Zrenner');
   catch
       save('tmp.mat')
   end
end
    
function [data,truePhase] = generateData(type, Fs,time, PAC, AAC)  

    switch type
        case 'sines_w_pink' 
            [pn] = make_pink_noise(1,time*Fs,1/Fs);
            Vlo = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % random movement in frequency, but smoothed so that it doesn't change faster than the oscillation rate
%             Vhi = (2+.25*randn(1,(time*Fs))) .* cos(2*pi*(50).*[1/Fs:1/Fs:time]);

            data = Vlo  + 10*pn;
            truePhase = angle(hilbert(Vlo))';
            
        case 'sines_w_white'
            Vlo = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % random movement in frequency, but smoothed so that it doesn't change faster than the oscillation rate
%             Vhi = (2+.25*randn(1,(time*Fs))) .* cos(2*pi*(50).*[1/Fs:1/Fs:time]);

            data = Vlo +  randn(1,length(Vlo));
            truePhase = angle(hilbert(Vlo))';

        case 'pink_noise'
            [pn] = make_pink_noise(1.5,time*Fs,1/Fs);
            AIC_comp = 8; % I think this is the number in the paper
            [XX,P,Vlo,Vhi,t] = simfun(0,0,'pink',0.05,'ci',AIC_comp);
            Vlo = 10*(Vlo/std(Vlo)); % forcing variance to 4
            data = (Vlo) + 10*pn; 
            truePhase = angle(hilbert(Vlo))';
            
        case 'SP_mdl'
            [y_origMdl, all_x_origMdl] = simSoulatModel([6], Fs,time,[.99],[10], 1);
            data =y_origMdl;
            truePhase = angle(all_x_origMdl(:,1) +1i*all_x_origMdl(:,2));
    end
end

function parsave(fname, structVar)
    save('-mat',fname, 'structVar')
end
