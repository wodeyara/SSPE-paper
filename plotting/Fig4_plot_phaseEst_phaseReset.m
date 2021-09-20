timePointsSplice = [3500, 4750, 6500,8500];
indsToTestPhase = [];
for i = 3
    indsToTestPhase = [indsToTestPhase, timePointsSplice(i)-167:timePointsSplice(i)+167];
end

% plotting Figure 4 of SSPE paper
figure
iter = 110;
plot(wrapToPi(simData(iter).truePhase(indsToTestPhase)), 'Linewidth',3)
hold on
plot(simData(iter).filt_phase(indsToTestPhase), 'Linewidth',1.25)
plot(simData(iter).AW_phase(indsToTestPhase), 'Linewidth',1.25)
plot(simData(iter).zrenner_phase(indsToTestPhase-999), 'Linewidth',1.25)
plot(simData(iter).SP_phase(indsToTestPhase), 'Linewidth',1.25)
legend({'True','acausal FIR','Blackwood','Zrenner','SSPE'})
% 
% plot([126,126], get(gca,'YLim'), 'Linewidth',1.25, 'color', 'r','Linestyle', '--')
% hold on
% plot([126,126], get(gca,'YLim'), 'Linewidth',1.25, 'color', 'r','Linestyle', '--')
% plot([251,251], get(gca,'YLim'), 'Linewidth',1.25, 'color', 'r','Linestyle', '--')
% plot([376,376], get(gca,'YLim'), 'Linewidth',1.25, 'color', 'r','Linestyle', '--')
% 
