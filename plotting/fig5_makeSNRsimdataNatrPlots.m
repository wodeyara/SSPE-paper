%% generate data
Fs = 1000;
timeAmt = 15;
L = Fs * timeAmt;
dt = 1/Fs;
alpha = 1;
ang_var2dev = @(v)(sqrt(-2*log(1-v)));

%%

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

cnt = 1;
for SNR = [1,2.5,5,7.5,10]
    x1 = randn(L,1);
    xf1 = fft(x1);
    A = abs(xf1);
    phase = angle(xf1);

    df = 1.0 / (dt*length(x1));
    faxis = (0:(length(x1)/2))*df;
    faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
    oneOverf = 1.0 ./ faxis.^alpha;
    oneOverf(1)=1;

    Anew = sqrt((A.^2).*oneOverf');
    xf1new = Anew .* exp(1i*phase);
    x1new = (ifft(round(xf1new,8)))';
    x1new = x1new * 1000;

    x1 = randn(L,1);
    xf1 = fft(x1);
    A = abs(xf1);
    phase = angle(xf1);

    df = 1.0 / (dt*length(x1));
    faxis = (0:(length(x1)/2))*df;
    faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
    oneOverf = 1.0 ./ faxis.^alpha;
    oneOverf(1)=1;

    Anew = sqrt((A.^2).*oneOverf');
    xf1new = Anew .* exp(1i*phase);
    x2new = (ifft(round(xf1new,8)))';
    x2new = x2new * 1000;

    lowAct = filtfilt(filtwts,1,x2new);
    % SNR as std dev ratio of signal/noise
    SNR_start = std(lowAct)/std(x1new);
    data(cnt,:) = (SNR/(SNR_start))*lowAct + x1new; % SNR = 1; sigmaFreqs [1,.05]
    data(cnt,:) = data(cnt,:)/std(data(cnt,:));
    truePhase(cnt,:) = angle(hilbert((SNR/(SNR_start))*lowAct));
    cnt = cnt + 1;
end

%%
cols = parula(5);
cols(5,:) = cols(5,:) - [.05,.05,.05];
SNR = [1,2.5,5,7.5,10];
for i = 1:5
    subplot(5,1,i)
    yyaxis left
    plot(1/Fs:1/Fs:3, data(i,1:3000), 'Linewidth',3,'Color',cols(i,:))
    set(gca,'Ylim', [-3,3])
    yyaxis right
    plot(1/Fs:1/Fs:3,truePhase(i,1:3000), 'Linewidth',2,'Color','blue')
    grid on
    title(['SNR = ',num2str(SNR(i))], 'Fontsize', 14)
    set(gca,'XTickLabels',{},'Fontsize',14)
end
%% making plots of the power spectrum

cols = parula(5);freqs = [1:50];

figure(1)
[pxx,f] = pmtm(data(1,:),5,freqs, 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(1,:))
hold on
[pxx,f] = pmtm(data(2,:),5,freqs, 1000)
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(2,:))
[pxx,f] = pmtm(data(3,:),5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(3,:))
[pxx,f] = pmtm(data(4,:),5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(4,:))
[pxx,f] = pmtm(data(5,:),5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(5,:))
set(gca,'XTick', [1,6,10,20,40,79], 'YTick', [-50:5:25], 'Ylim',[-50,25])

grid
ylabel('Power (dB)')
xlabel('Frequency (Hz)')
set(gca,'Fontsize',13)
legend({'SNR = 1','SNR = 2.5','SNR = 5','SNR = 7.5','SNR = 10'})
ylim([-50,0])
