function [init_f,init_a,init_sigma,R0] = initializeParams(y0,num_component,T1,T2)
    % borrowed code from the MK model approach to initialize parameters
    % appropriately.
    % Last Edit: Ani Wodeyar 2/1/2020

    T0 = length(y0);
    y = y0(T1:T2);
    T = length(y);
    MAX_AR = 2*num_component+1;
    AR_fit;

    minAIC = inf;
    minAIC2 = inf;
    for ARdeg=num_component:2*num_component
        tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
        if ARdeg-nnz(imag(tmp))/2 == num_component && AIC_ARwithnoise(ARdeg) < minAIC
            z0 = tmp;
            E0 = ARwithnoise_param(ARdeg,ARdeg+1);
            R0 = ARwithnoise_param(ARdeg,ARdeg+2);
            minAIC = AIC_ARwithnoise(ARdeg);
            optARdeg = ARdeg;
        end
        if ARdeg-nnz(imag(tmp))/2 >= num_component && AIC_ARwithnoise(ARdeg) < minAIC2
            z1 = tmp;
            E1 = ARwithnoise_param(ARdeg,ARdeg+1);
            R1 = ARwithnoise_param(ARdeg,ARdeg+2);
            minAIC2 = AIC_ARwithnoise(ARdeg);
            optARdeg2 = ARdeg;
        end
    end
    if minAIC == inf
        if minAIC2 == inf
            warning('no AR model with %d oscillators',num_component);
        end
        z0 = z1;
        E0 = E1;
        R0 = R1;
        optARdeg = optARdeg2;
        [B,I] = sort(abs(z0),'descend');
        z0 = z0(I);
    end
    [tmp,ARdeg] = min(AIC_ARwithnoise);
    tmp = roots([1 ARwithnoise_param(ARdeg,1:ARdeg)]);
    if nnz(imag(tmp)>=0) >= num_component
        z0 = tmp;
        E0 = ARwithnoise_param(ARdeg,ARdeg+1);
        R0 = ARwithnoise_param(ARdeg,ARdeg+2);
        optARdeg = ARdeg;
    end
    VV = zeros(optARdeg,optARdeg);
    for j=1:optARdeg
        for i=1:optARdeg
            VV(i,j) = z0(j)^(1-i);
        end
    end
    QQ = inv(VV)*[E0 zeros(1,optARdeg-1); zeros(optARdeg-1,optARdeg)]*inv(VV)';
    [B,I] = sort(diag(real(QQ))./(1-abs(z0).^2),'descend');
    z0 = z0(I);

%     init_a = [];
%     init_f = [];
%     z_complex = (imag(z0) ~=0);
%     if sum(z_complex) >0
%         init_a = (abs(z0(z_complex)));
%         init_f = unique((1/(2*pi))*acos(abs(real(z0(z_complex)))./init_a));
%         init_a = unique(init_a);
%     end
%     % covering cases where we have no imaginary root
%     if length(init_a) < num_component
%        init_f = [init_f, zeros(1, num_component - length(init_f))];
%        tmp  =abs(z0(~z_complex));
%        init_a = [init_a, tmp(1:num_component - length(init_a))];
%     end
    
    init_a = zeros(1,num_component);
    init_f = zeros(1,num_component);
    kk = 1;
    for k=1:num_component
        init_a(k) = abs(z0(kk));
        init_f(k) = abs(angle(z0(kk)));
        if isreal(z0(kk)) && z0(kk) < 0
            init_f(k) = 4; % for pi
        end
        if imag(z0(kk)) == 0
            kk = kk+1;
        else
            kk = kk+2;
        end
    end

    [B,I] = sort(init_f);
    init_a = init_a(I);
    init_f = init_f(I);

    if mod(T,2) == 0
        freq = [2*pi/T*(0:T/2-1) pi];
    else
        freq = [2*pi/T*(0:(T-1)/2)];
    end
    init_f(init_f==4) = pi;
    P = zeros(length(freq),num_component);
    for k=1:num_component
        a = init_a(k);
        theta = init_f(k);
        A = (1-2*a^2*cos(theta)^2+a^4*cos(2*theta))/a/(a^2-1)/cos(theta);
        b = (A-2*a*cos(theta)+sign(cos(theta))*sqrt((A-2*a*cos(theta))^2-4))/2;
        for j=1:length(freq)
            P(j,k) = -a*cos(theta)/b*abs(1+b*exp(-1i*freq(j)))^2/abs(1-2*a*cos(theta)*exp(-1i*freq(j))+a^2*exp(-2*1i*freq(j))).^2;
        end
    end
    p = zeros(length(freq),1);
    for j=1:length(freq)
        p(j) = abs(y*exp(-1i*freq(j)*(0:T-1)'))^2/T;
    end
    if cond(P'*P) < 10^6    % modified 2017/01/14
        init_sigma = (P'*P)\(P'*p);
    else
        init_sigma = expGLMfit(P,p);
    end
    init_sigma(init_sigma<0) = R0;