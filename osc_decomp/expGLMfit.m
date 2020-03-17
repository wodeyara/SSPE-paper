function beta = expGLMfit(X,y)
    options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000);
    beta = exp(fminunc(@(b)expGLM_ll(X,y,exp(b)),zeros(size(X,2),1),options));
end

