%% Figure 6
subplot(212)
params.tapers = [10,19];
params.Fs = 1000;
params.pad = -1;
params.fpass = [1,30];

[S,T, F]  = mtspecgramc([zeros(4000,1);data(5e4+1:30e4);zeros(5000,1)],[10,1],params);
subplot(212)
T = T(12:246)-5.001;
spect = 10*log10(abs(S(12:246,:))');
imagesc(T,F,spect)
set(gca,'ydir','norm')
caxis([15,35])

xlabel('Time (s)')
set(gca,'Fontsize', 16)
ylabel('Frequency')
h = colorbar;
ylabel(h,'Power (dB)')

subplot(211)
errorEst(errorEst == 0) = NaN;
fifthAndNinetyFifth = prctile(errorEst(:,51:end)',[5,95]);
tmp = fifthAndNinetyFifth;
tmp2(1:2,1) = tmp(1:2,1);
for i = 2:size(tmp,2)-1
    if isnan(tmp(1,i))&& ~isnan(tmp(1,i-1))
       tmp2(1:2,i) = mean(tmp(1:2,i-1));
    elseif isnan(tmp(1,i))&& ~isnan(tmp(1,i+1))
        tmp2(1:2,i) = mean(tmp(1:2,i+1));
    else
        tmp2(1:2,i) = tmp(1:2,i);
    end
end
tmp2(1:2,size(tmp,2)) = tmp(1:2,end);
% tmp(isnan(tmp)) = 0;
shade(11:245, tmp2(1,11:245), 11:245,tmp2(2,11:245) , 'FillType',[2 1], 'Color', 'red', 'FillAlpha', .65, 'FillColor', [.9,.1,.1])
hold on
scatter(11:245,errorEst(11:245,51), 'filled', 'blue')
plot(11:245,errorEst(11:245,51),'Linewidth',2,'Color','blue')
xlim([11,245])
set(gca,'Fontsize', 13)
xlabel('Time (sec)')
ylabel('Circular Deviation (deg)')

%% Figure 8

subplot(212)
params.tapers = [10,19];
params.Fs = 1000;
params.pad = -1;
params.fpass = [1,30];
[S,T, F]  = mtspecgramc([zeros(4000,1);data(5e4+1:25e4)';zeros(5000,1)],[10,1],params);
T = T(12:196)-5.001;
spect = 10*log10(abs(S(12:196,:))');
imagesc(T,F,spect)
set(gca,'ydir','norm')
caxis([-15,-.5])
xlabel('Time (s)')
set(gca,'Fontsize', 13)
ylabel('Frequency')
h = colorbar;
ylabel(h,'Power (dB)')
xlim([11,195])


subplot(211)
errorEst(errorEst == 0) = NaN;

fifthAndNinetyFifth = prctile(errorEst(:,51:end)',[2.5,97.5]);
tmp = fifthAndNinetyFifth;
tmp(isnan(tmp)) = 0;
shade(11:195, tmp(1,11:195), 11:195,tmp(2,11:195) , 'FillType',[2 1], 'Color', 'red', 'FillAlpha', .65, 'FillColor', [.9,.1,.1])
hold on
scatter(11:195,errorEst(11:195,51), 'filled', 'blue')
plot(11:195,errorEst(11:195,51),'Linewidth',2,'Color','blue')
xlim([11,195])
set(gca,'Fontsize', 13)
xlabel('Time (sec)')
ylabel('Circular Deviation (deg)')

%% getting PSD
tmp1 = logical([zeros(4,1);ones(191,1);zeros(5,1)]);
allColors = parula(10);

% S = 10*log10(S);
allPSD_mean = mean(S(tmp1,:),1);
allPSD_std = (1/sqrt(size(S(tmp1,:),1))) * std(S);

h = plot(F,10*log10(allPSD_mean), 'linewidth', 2, 'Color', allColors(1,:));
hold on
g = shade(F, 10*log10(allPSD_mean - 1.96*allPSD_std)...
    , F, 10*log10(allPSD_mean + 1.96*allPSD_std) ,...
    'FillType',[2 1], 'FillAlpha', .65, 'FillColor', allColors(3,:), 'Color', allColors(3,:));
set(get(get(g(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(g(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(g(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'Fontsize', 16)
ylabel('Power (dB)')
xlabel('Frequency')
grid on
xlim([1,50])
