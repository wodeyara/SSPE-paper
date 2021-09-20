options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000 ,'UseParallel', ~0);
AIC_ARraw = zeros(1,MAX_AR);
ARraw_param = zeros(MAX_AR,MAX_AR+1);
AIC_ARwithnoise = zeros(1,MAX_AR);
ARwithnoise_param = zeros(MAX_AR,MAX_AR+2);
for ARdeg=1:MAX_AR
    [A,E] = AR_MLE(y,ARdeg);
    a = zeros(ARdeg,ARdeg);
    a(:,ARdeg) = -A(2:ARdeg+1);
    for m=ARdeg:-1:2
        for i=1:m-1
            a(i,m-1) = (a(i,m)+a(m,m)*a(m-i,m))/(1-a(m,m)^2);
        end
    end
    c = zeros(ARdeg,1);
    for m=1:ARdeg
        c(m) = a(m,m);
    end
    AIC_ARraw(ARdeg) = 2*AR_ll(y,log(1+c)-log(1-c))+2*(ARdeg+1);
    ARraw_param(ARdeg,1:ARdeg+1) = [A(2:ARdeg+1) E];
    z = roots(A);
    k = 1;
    while max(abs(z)) >= 1
        [A,E] = AR_MLE(y,ARdeg-k);
        z = roots(A);
        A = [A zeros(1,k)];
        k = k+1;
    end
    [A,E,R] = armyule(y,ARdeg);
    param = [A(2:ARdeg+1) log(E)-log(R)];
    [mll,R] = ARwithnoise_ll(y,param);
    AIC_ARwithnoise(ARdeg) = 2*mll+2*(ARdeg+2);
    ARwithnoise_param(ARdeg,1:ARdeg+2) = [param(1:ARdeg) exp(param(ARdeg+1))*R R];
end
