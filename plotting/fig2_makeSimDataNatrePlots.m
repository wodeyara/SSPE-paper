cols = jet(24);
cols = cols(1:6:24,:);

time = 10;
Fs = 1000;
[pn] = make_pink_noise(1.5,time*Fs,1/Fs);
Vlo = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % 
truePhase = wrapToPi(2*pi*(6).*[1/Fs:1/Fs:time]);
data = Vlo + 10*pn;
subplot(412)
% plot(1/Fs:1/Fs:2,data(1:2000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
% hold on
yyaxis left
% plot(1/Fs:1/Fs:2, Vlo(1:2000),'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',2)
% hold on
plot(1/Fs:1/Fs:2, data(1:2000), 'Linewidth',2,'Color',cols(2,:))
yyaxis right
plot(1/Fs:1/Fs:2,truePhase(1:2000), 'Linewidth',2,'Color','blue')

title('Sines In Pink Noise','Fontsize', 16)
set(gca,'XTickLabels',{},'Fontsize',14)

Vlo = (10).*cos(2*pi*(6).*[1/Fs:1/Fs:time]); % 
data = Vlo + randn(1,length(Vlo));
truePhase = wrapToPi(2*pi*(6).*[1/Fs:1/Fs:time]);
subplot(411)
% plot(1/Fs:1/Fs:2,data(1:2000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
% hold on
yyaxis left
% plot(1/Fs:1/Fs:2, Vlo(1:2000),'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',2)
% hold on
plot(1/Fs:1/Fs:2, data(1:2000), 'Linewidth',2,'Color',cols(1,:))
yyaxis right
plot(1/Fs:1/Fs:2,truePhase(1:2000), 'Linewidth',2,'Color','blue')
title('Sines in White Noise', 'Fontsize', 16)
set(gca,'XTickLabels',{},'Fontsize',14)

[pn] = make_pink_noise(1.5,10*Fs,1/Fs);
AIC_comp = 8; % I think this is the number in the paper
[XX,P,Vlo,Vhi,t] = simfun(0,0,'pink',0.05,'ci',AIC_comp);
Vlo = 10*(Vlo/std(Vlo)); % forcing variance to 4
data = (Vlo) + 10*pn(1:1e4); 
truePhase = angle(hilbert(Vlo))';
subplot(413)
% plot(1/Fs:1/Fs:2,data(1:2000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
% hold on
yyaxis left
% plot(1/Fs:1/Fs:2, Vlo(1:2000),'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',2)
% hold on
plot(1/Fs:1/Fs:2, data(1:2000), 'Linewidth',2,'Color',cols(3,:))
yyaxis right
plot(1/Fs:1/Fs:2,truePhase(1:2000), 'Linewidth',2,'Color','blue')

title('Filtered Pink Noise','Fontsize', 16)
set(gca,'XTickLabels',{},'Fontsize',14)

[y_origMdl, all_x_origMdl] = simSoulatModel([6], Fs,time,[.99],[10], 1);
data =y_origMdl;
truePhase = angle(all_x_origMdl(:,1) +1i*all_x_origMdl(:,2));
subplot(414)
% plot(1/Fs:1/Fs:2,data(1:2000),'Color',[0.8500, 0.3250, 0.0980],'Linewidth',2)
% hold on
yyaxis left
% plot(1/Fs:1/Fs:2, all_x_origMdl(1:2000,1),'Color',[0.4940, 0.1840, 0.5560], 'Linewidth',2)
% hold on
plot(1/Fs:1/Fs:2, data(1:2000), 'Linewidth',2,'Color',cols(4,:))
% legend({'6 Hz Observation','Phase'})
yyaxis right
plot(1/Fs:1/Fs:2,truePhase(1:2000), 'Linewidth',2,'Color','blue')

title('MK model','Fontsize', 16)
set(gca,'Fontsize',14)
