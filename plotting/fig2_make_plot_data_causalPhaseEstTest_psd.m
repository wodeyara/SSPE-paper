freqs = [1:50];

cols = jet(24);
cols = cols(1:6:24,:);

[pn] = make_pink_noise(1.5,time*Fs,1/Fs);
Vlo1 = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % 
data2 = Vlo1 + 10*pn;
% plot(1/Fs:1/Fs:2,data(1:2000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
% hold on

Vlo2 = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % 
data1 = Vlo2 +  randn(1,length(Vlo1));

[pn] = make_pink_noise(1.5,10*Fs,1/Fs);
AIC_comp = 8; % I think this is the number in the paper
[XX,P,Vlo,Vhi,t] = simfun(0,0,'pink',0.05,'ci',AIC_comp);
Vlo3 = 10*(Vlo/std(Vlo)); % forcing variance to 4
data3 = (Vlo3) + 10*pn(1:1e4); 

[y_origMdl, all_x_origMdl] = simSoulatModel([6], Fs,time,[.99],10, 1);
data4 =y_origMdl;
Vlo4 = all_x_origMdl(:,1);
truePhase = angle(all_x_origMdl(:,1) +1i*all_x_origMdl(:,2));

figure(1)
subplot(1,2,1)
[pxx,f] = pmtm(Vlo1,5,freqs, 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(1,:))
hold on
[pxx,f] = pmtm(Vlo2,5,freqs, 1000)
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(2,:))
[pxx,f] = pmtm(Vlo3,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(3,:))
[pxx,f] = pmtm(Vlo4,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(4,:))
set(gca,'XTick', [1,6,10,20,40,79], 'YTick', [-50:5:25], 'Ylim',[-50,25])

grid
ylabel('Power/frequency (dB/Hz)')
xlabel('Frequency (Hz)')
title('Signal')
set(gca,'Fontsize',13)
legend({'Sines with White Noise', 'Sines in Pink Noise', 'Filtered Pink Noise','MK model'})

subplot(1,2,2)
[pxx,f] = pmtm(data1,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(1,:))
hold on
[pxx,f] = pmtm(data2,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(2,:))
[pxx,f] = pmtm(data3,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(3,:))
[pxx,f] = pmtm(data4,5,[freqs], 1000);
plot(f,10*log10(pxx),'Linewidth',1.5, 'Color', cols(4,:))
set(gca,'XTick', [1,6,10,20,40,79], 'YTick', [-50:5:25], 'Ylim',[-50,25])

grid

ylabel('Power/frequency (dB/Hz)')
xlabel('Frequency (Hz)')
title('Observation')
set(gca,'Fontsize',13)
legend({'Sines with White Noise', 'Sines in Pink Noise', 'Filtered Pink Noise','MK model'})
