function []=plot_phase(y,fs,phi_prop,decomp_mu,decomp_cov,num_component)
    nmc = 10^4;
    T = length(y);
    T1 = 1;
    T2 = T;
    Tstart = 0;
    phase1 = zeros(1,T2-T1+1);
    phase2 = zeros(1,T2-T1+1);
    seeds = randn(2,nmc);
    figure;
    for k=1:num_component
        for t=T1:T2
            tmp = angle_conf_MC(decomp_mu(2*k-1:2*k,t,num_component),decomp_cov(2*k-1:2*k,2*k-1:2*k,t,num_component),2*normcdf(1)-1,seeds);
            phase1(t-T1+1) = tmp(1);
            phase2(t-T1+1) = tmp(2);
        end
        subplot(num_component,1,k)
        hold on
        plot_phase_area(Tstart+(T1:T2)/fs,phi_prop(k,T1:T2,num_component),phase1,phase2,[.8,.8,.8]);
        plot_phase_nocross(Tstart+(T1:T2)/fs,phi_prop(k,T1:T2,num_component),'b',2);
        xlim(Tstart+[T1 T2]/fs);
        set(gca,'FontSize',12);
        ax = gca;
        set(ax,'YTick',[-pi 0 pi]);
        set(ax,'YTickLabel',{'-3.14', '0', '3.14'});
    end
end
