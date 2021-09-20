function [phi_prop,decomp_mu,decomp_cov,AIC_osc,osc_param] = osc_decomp_fast(y0,fs,MAX_COMPONENT,T1,T2)
%
% fast version of osc_decomp
%
% This function reduces computational cost by using only part of the input time series for parameter estimation
%
% Input:
%    y0:            time series (1 times T0)
%    fs:            sampling frequency (scalar)
%    MAX_COMPONENT: maximum number of oscillation components (scalar)
%    T1:            first time point of the interval for parameter estimation (scalar)
%    T2:            last time point of the interval for parameter estimation (scalar)
%                   (Namely, y0(T1:T2) is used for parameter estimation)
%
% Output:
%    phi_prop:      estimated phase of each oscillation component (MAX_COMPONENT times T0 times MAX_COMPONENT)
%                   phi_prop(k,:,K) is the estimated phase of the k-th
%                   oscillator in the decomposition into K components
%    decomp_mu:     smoothed coordinate of each oscillator (2*MAX_COMPONENT times T0 times MAX_COMPONENT)
%                   decomp_mu(2*k-1:2*k,:,K) is the smoothed coordinate of the k-th oscillator in the decomposition into K components
%    decomp_cov:    smoothed covariance of each oscillator (2*MAX_COMPONENT times 2*MAX_COMPONENT times T0 times MAX_COMPONENT)
%                   decomp_cov(2*k-1:2*k,2*k-1:2*k,:,K) is the smoothed covariance of the k-th oscillator in the decomposition into K components
%    AIC_osc:       AIC of the oscillator model (1 times MAX_COMPONENT)
%                   AIC_osc(K) is the AIC of the oscillator model with K oscillation components
%    osc_param:     estimated parameters of the oscillator model (MAX_COMPONENT times 3*MAX_COMPONENT+1)
%                   osc_param(K,1:K) is the estimated a_1,...,a_K of the oscillator model with K oscillation components
%                   osc_param(K,K+1:2*K) is the estimated f_1,...f_K of the oscillator model with K oscillation components
%                   osc_param(K,2*K+1:3*K) is the estimated sigma_1^2,...,sigma_K^2 of the oscillator model with K oscillation components
%                   osc_param(K,3*K+1) is the estimated tau^2 of the oscillator model with K oscillation components
%
    options = optimoptions('fminunc','Algorithm','quasi-newton','MaxFunEvals',10000,'UseParallel', ~0 );
    T0 = length(y0);
    phi_prop = zeros(MAX_COMPONENT,T0,MAX_COMPONENT);
    decomp_mu = zeros(2*MAX_COMPONENT,T0,MAX_COMPONENT);
    decomp_cov = zeros(2*MAX_COMPONENT,2*MAX_COMPONENT,T0,MAX_COMPONENT);
    AIC_osc = zeros(1,MAX_COMPONENT);
    osc_param = zeros(MAX_COMPONENT,3*MAX_COMPONENT+1);
    y = y0(T1:T2);
    T = length(y);
    MAX_AR = 2*MAX_COMPONENT;
    AR_fit;
    for num_component=1:MAX_COMPONENT
        num_component
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
%    	nf0 = nnz(init_f==0);
%    	nf1 = nnz(init_f==4);
%	    freq = init_f;
%    	freq(2:nf0) = pi*rand(1,nf0-1);
%	    freq(num_component-nf1+2:num_component) = pi*rand(1,nf1-1);
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
        
        param = fminunc(@(param)prop_ll(y,param,init_f),[atanh(2*init_a-1) zeros(1,num_component) log(init_sigma'/R0)],options);
        [mll,sigman] = prop_ll(y,param,init_f);
        param(num_component+1:2*num_component) = abs(init_f+tanh(param(num_component+1:2*num_component))*pi);
        AIC_osc(num_component) = 2*mll+2*(3*num_component+1);
	    [tmp,I] = sort(param(num_component+1:2*num_component));
        a = (tanh(param(I))+1)/2;
        theta = param(num_component+I);
        sigma = exp(param(2*num_component+I))*sigman;
        osc_param(num_component,1:3*num_component+1) = [a theta*fs/2/pi sigma sigman];
    
        m = 2*num_component;
        x_pred1 = zeros(m,T0);
        x_filt = zeros(m,T0);
        x_smooth = zeros(m,T0);
        V_pred1 = zeros(m,m,T0);
        V_filt = zeros(m,m,T0);
        V_smooth = zeros(m,m,T0);
        F = zeros(m,m);
        Q = zeros(m,m);
        H = zeros(1,m);
        for k=1:num_component
            F(2*k-1:2*k,2*k-1:2*k) = a(k)*[cos(theta(k)) -sin(theta(k)); sin(theta(k)) cos(theta(k))];
            Q(2*k-1,2*k-1) = 1-a(k)^2;
            Q(2*k,2*k) = 1-a(k)^2;
        end
        H(1:2:m) = sqrt(sigma)./sqrt(1-a.^2);
        R = sigman;
        x_pred1(:,1) = zeros(m,1);
        for k=1:num_component
            V_pred1(2*k-1:2*k,2*k-1:2*k,1) = eye(2);
        end
        for t=1:T0-1
            Kg = V_pred1(:,:,t)*H'*inv(H*V_pred1(:,:,t)*H'+R);
            x_filt(:,t) = x_pred1(:,t) + Kg*(y0(t)-H*x_pred1(:,t));
            V_filt(:,:,t) = (eye(m)-Kg*H)*V_pred1(:,:,t);
            x_pred1(:,t+1) = F*x_filt(:,t);
            V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
        end
        Kg = V_pred1(:,:,T0)*H'*inv(H*V_pred1(:,:,T0)*H'+R);
        x_filt(:,T0) = x_pred1(:,T0) + Kg*(y0(T)-H*x_pred1(:,T0));
        V_filt(:,:,T0) = (eye(m)-Kg*H)*V_pred1(:,:,T0);
        x_smooth(:,T0) = x_filt(:,T0);
        V_smooth(:,:,T0) = V_filt(:,:,T0);
        for t=T0-1:-1:1
            x_smooth(:,t) = x_filt(:,t) + V_filt(:,:,t)*F'*(V_pred1(:,:,t+1)\(x_smooth(:,t+1)-x_pred1(:,t+1)));
            CC = V_filt(:,:,t)*F'*inv(V_pred1(:,:,t+1));
            V_smooth(:,:,t) = V_filt(:,:,t) + CC*(V_smooth(:,:,t+1)-V_pred1(:,:,t+1))*CC';
        end
        for k=1:num_component
            phi_prop(k,:,num_component) = atan2(x_smooth(2*k,:),x_smooth(2*k-1,:));
        end
        decomp_mu(1:m,:,num_component) = x_smooth;
        decomp_cov(1:m,1:m,:,num_component) = V_smooth;
    end
end


