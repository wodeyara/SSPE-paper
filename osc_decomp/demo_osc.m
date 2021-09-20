load('WolferSunspotData');
y = log(sunspot+1);
fs = 1;
MAX_COMPONENT = 5;
tic
[decomp_phase,decomp_mu,decomp_cov,AIC_osc,osc_param] = osc_decomp(y,fs,MAX_COMPONENT);
time2=toc
[tmp,K] = min(AIC_osc);
plot_decomp(y,fs,decomp_mu,decomp_cov,osc_param,K);
plot_phase(y,fs,decomp_phase,decomp_mu,decomp_cov,K);
osc_a = osc_param(K,1:K);
osc_f = osc_param(K,K+1:2*K);
osc_sigma = osc_param(K,2*K+1:3*K);
osc_r = osc_param(K,3*K+1);
str = sprintf('The number of oscillators is K=%d.',K);
disp(str);
str = 'The periods of K oscillators are ';
for k=1:K-1
    str = [str sprintf('%.2f, ',1./osc_f(k))];
end
str = [str sprintf('%.2f years.',1./osc_f(K))];
disp(str);
