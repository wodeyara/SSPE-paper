

tmpPhase = phase(250,:);
tmpPhase = round(10*tmpPhase/5)*5/10;

subplot(311)
violinplot(log10(mep_amplitude),tmpPhase,parula(20))
grid on
set(gca,'Fontsize', 13)
SSPE_mdl = fitlm([sin(phase(250,:))', cos(phase(250,:))'],log10(mep_amplitude));
Ypred = predict(SSPE_mdl, [sin([-pi+pi/100:pi/100:pi])',cos([-pi+pi/100:pi/100:pi])']);
plot(.5+ (13/length(Ypred)):13/length(Ypred):13.5 , Ypred,'Linewidth',2)
title(['No CI threshold: ','R^2 = ', num2str(round(SSPE_mdl.Rsquared.Ordinary,2))]);
inds = find((phaseWidth(250,:) < 50));

subplot(312)
violinplot(log10(mep_amplitude(inds)),tmpPhase(inds),parula(20))
grid on
set(gca,'Fontsize', 13)
SSPE_mdl = fitlm([sin(phase(250,inds))', cos(phase(250,inds))'],log10(mep_amplitude(inds)));
title(['CI<50 degrees: ','R^2 = ', num2str(round(SSPE_mdl.Rsquared.Ordinary,2))]);
Ypred = predict(SSPE_mdl, [sin([-pi+pi/100:pi/100:pi])',cos([-pi+pi/100:pi/100:pi])']);
plot(.5+ (13/length(Ypred)):13/length(Ypred):13.5 , Ypred,'Linewidth',2)

inds = find((phaseWidth(250,:) < 25));
subplot(313)
violinplot(log10(mep_amplitude(inds)),tmpPhase(inds),parula(20))
SSPE_mdl = fitlm([sin(phase(250,inds))', cos(phase(250,inds))'],log10(mep_amplitude(inds)));
hold on
Ypred = predict(SSPE_mdl, [sin([-pi+pi/100:pi/100:pi])',cos([-pi+pi/100:pi/100:pi])']);
plot(.5+ (13/length(Ypred)):13/length(Ypred):13.5 , Ypred,'Linewidth',2)
title(['CI<25 degrees: ','R^2 = ', num2str(round(SSPE_mdl.Rsquared.Ordinary,2))]);
grid on
set(gca,'Fontsize', 13)
xlabel('Phase')
ylabel('log_{10}(MEP amplitude)')