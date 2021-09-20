%%
% this was set up for a local field potential with theta band. The code
% below plots the amplitude against the credible interval - see Figure 7C.

%%
[map] = brewermap(6,'Dark2');
subplot(211)
% h = scatter(log10(credInt(2:end)), log10(amp(2:end)), 50, [.5,.5,.5],'filled','HandleVisibility','off');
% h.MarkerFaceAlpha = .2;
A = log10(credInt);
B = log10(amp);
AB =[A B];
% Find unique rows and corresponding indices
[uniqueAB, ~, n] = unique(round(AB,2), 'rows');
% Find number of occurrences
nHist = hist(n, unique(n));
nHist = log10(nHist);
mx = max(floor(100*(nHist + abs(min(nHist)))))+1;
% Create colors for each number of occurrence
colors = jet(mx);
colormap(colors);
% Construct a color matrix
cMatrix = colors(floor(100*(nHist + abs(min(nHist))))+1, :);
% Create scatter plot
scatter(uniqueAB(:, 1), uniqueAB(:, 2), 25, cMatrix, 'filled');
% colorbar
colorbar('YTick', [.01, .5,.99], ...
         'YTickLabel', round(10.^(round(prctile(nHist,[0,50,99])))));
hold on

cnt = 1;
% ampThrEnd = [0,15,45,75,100,150]; 
for ampThr = prctile(amp,[95,80,65]) %[15,45,75,100,150]
    
%     yyaxis right
%     [p,X]= histcounts(credInt(amp > ampThr),'normalization','pdf');
%     plot(X(2:end),p/sum(p),'linewidth',2, 'Color', map(cnt,:))
%     xlim([0,pi])
%     hold on
    
%     yyaxis left
    plot(get(gca,'Xlim'), [log10(ampThr),log10(ampThr)],...
        'linewidth',3,'Color',map(cnt,:));
    cnt = cnt + 1;
end
% xlim([-3,log10(pi)])
xlim([-1.8,.5])
legend({'95th percentile','80th percentile','65th percentile'})
set(gca,'Fontsize', 16)
ylabel('log_{10}(Amplitude \muV)' )

subplot(212)
credIntRed = ([credInt(amp>prctile(amp,95))',credInt(amp>prctile(amp,80))', ...
    credInt(amp>prctile(amp,65))']);
cats = [ones(1,nnz(amp>prctile(amp,95))),2*ones(1,nnz(amp>prctile(amp,80))),...
        3*ones(1,nnz(amp>prctile(amp,65)))];
violinplot(credIntRed,cats,map(1:3,1:3),'ViolinAlpha',.05) 
set(gca,'YScale','log')
set(gca,'XTickLabel',{'95th','80th','65th'})
set(gca,'Fontsize', 16)
ylabel('Credible Interval')
ylim([10^(-1.8),10^(.5)])
view([90,90])