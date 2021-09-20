% Plotting Figure 2 of SSPE paper
%%
ang_var2dev = @(v) sqrt(-2*log(1-v));
y = zeros(1000, 4,4);
load('causalPhaseEst_sim_sines_w_white.mat', 'var_FIR_Hilb', ...
                    'var_AW','var_SPcausal','var_Zrenner')
y(:,1,1) = rad2deg(ang_var2dev(var_FIR_Hilb));
y(:,1,2) = rad2deg(ang_var2dev(var_AW));
y(:,1,3) = rad2deg(ang_var2dev(var_Zrenner));
y(:,1,4) = rad2deg(ang_var2dev(var_SPcausal));

load('causalPhaseEst_sim_SPmdl.mat', 'var_FIR_Hilb', ...
                    'var_AW','var_SPcausal','var_Zrenner')
y(:,4,1) = rad2deg(ang_var2dev(var_FIR_Hilb));
y(:,4,2) = rad2deg(ang_var2dev(var_AW));
y(:,4,3) = rad2deg(ang_var2dev(var_Zrenner));
y(:,4,4) = rad2deg(ang_var2dev(var_SPcausal));

load('causalPhaseEst_sim_sines_w_pink.mat', 'var_FIR_Hilb', ...
                    'var_AW','var_SPcausal','var_Zrenner')
y(:,2,1) = rad2deg(ang_var2dev(var_FIR_Hilb));
y(:,2,2) = rad2deg(ang_var2dev(var_AW));
y(:,2,3) = rad2deg(ang_var2dev(var_Zrenner));
y(:,2,4) = rad2deg(ang_var2dev(var_SPcausal));

load('causalPhaseEst_sim_pink_noise.mat', 'var_FIR_Hilb', ...
                    'var_AW','var_SPcausal','var_Zrenner')
y(:,3,1) = rad2deg(ang_var2dev(var_FIR_Hilb));
y(:,3,2) = rad2deg(ang_var2dev(var_AW));
y(:,3,3) = rad2deg(ang_var2dev(var_Zrenner));
y(:,3,4) = rad2deg(ang_var2dev(var_SPcausal));

tmp = parula(4);
tmpCol1 = [.65,.65,.65];
tmpCol2 = tmp(2,:);
tmpCol3 = tmp(3,:);
tmpCol4 = tmp(4,:);
x = [1:4];
figure;
h = iosr.statistics.boxPlot(x,y,...
'symbolColor','k',...
'medianColor','k',...
'symbolMarker',{'+','o','d','x'},...
'boxcolor',{tmpCol1;tmpCol2;tmpCol3;tmpCol4},...
'groupLabels',{'acausal FIR','Blackwood et al. 2018','Zrenner et al. 2020','SSPE'},...
'showLegend',true);
box on
set(gca,'Fontsize', 16)
set(gca,'XTickLabel', {'Sines In White Noise)', 'Sines in Pink Noise', 'Filtered Pink Noise', 'State Space Model'})
ylabel('Circular Standard Deviation')
% ylim([0,.08])
