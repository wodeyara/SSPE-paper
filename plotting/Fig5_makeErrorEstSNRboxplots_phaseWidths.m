y = zeros(1000, 5,4);

y(:,:,1) = circVarFIR';
y(:,:,2) = circVarAWMeth';
y(:,:,3) = circVarZrenner';
y(:,:,4) = circVarMKMeth';

tmp = parula(4);
tmpCol1 = [.65,.65,.65];
tmpCol2 = tmp(2,:);
tmpCol3 = tmp(3,:);
tmpCol4 = tmp(4,:);
x = [1:5];
figure;
h = iosr.statistics.boxPlot(x,y,...
'symbolColor','k',...
'medianColor','k',...
'symbolMarker',{'+','o','d','x'},...
'boxcolor',{tmpCol1;tmpCol2;tmpCol3;tmpCol4},...
'groupLabels',{'acausal FIR','Blackwood','Zrenner','SSPE'},...
'showLegend',true);
box on
set(gca,'Fontsize', 16)
set(gca,'XTickLabel', {'SNR = 1', 'SNR = 2.5', 'SNR = 5', 'SNR = 7.5', 'SNR = 10'})
ylabel('Circular Standard Deviation')
% ylim([0,.08])

%%
% figure
% subplot(121)
violinplot(circStdDev,{'1', '2.5', '5', '7.5'}, parula(4))
xlabel('SNR')
view([90,90])
set(gca,'Fontsize',16)
% ylabel('Amplitude')
ylabel('Credible Interval Width')
grid on
% 
% subplot(122)
% violinplot(amps,{'1', '2.5', '5', '7.5', '10'}, parula(5))
% xlabel('SNR')
% view([90,90])
% set(gca,'Fontsize',16)
% ylabel('Amplitude')
% grid on