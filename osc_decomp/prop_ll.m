function [mll,Rhat] = prop_ll(y,param,init_f)
    T = length(y);
    num_osc = length(param)/3;
    param(num_osc+1:2*num_osc) = init_f+tanh(param(num_osc+1:2*num_osc))*pi;
    F = zeros(2*num_osc,2*num_osc);
    Q = zeros(2*num_osc,2*num_osc);
    for i=1:num_osc
        F(2*i-1:2*i,2*i-1:2*i) = (tanh(param(i))+1)/2*[cos(param(num_osc+i)) -sin(param(num_osc+i)); sin(param(num_osc+i)) cos(param(num_osc+i))];
        Q(2*i-1:2*i,2*i-1:2*i) = exp(param(2*num_osc+i))*eye(2);
    end
    H = zeros(1,2*num_osc);
    H(1:2:2*num_osc) = 1;
    R = 1;

    x_pred1 = zeros(2*num_osc,T);
    x_filt = zeros(2*num_osc,T);
    V_pred1 = zeros(2*num_osc,2*num_osc,T);
    V_filt = zeros(2*num_osc,2*num_osc,T);
    x_pred1(:,1) = zeros(2*num_osc,1);
    for i=1:num_osc
        V_pred1(2*i-1:2*i,2*i-1:2*i,1) = Q(2*i-1:2*i,2*i-1:2*i)/(1-F(2*i-1,2*i-1)^2-F(2*i-1,2*i)^2);
    end
    for t=1:T-1
        x_filt(:,t) = x_pred1(:,t) + V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\(y(t)-H*x_pred1(:,t)));
        V_filt(:,:,t) = V_pred1(:,:,t) - V_pred1(:,:,t)*H'*((H*V_pred1(:,:,t)*H'+R)\H)*V_pred1(:,:,t);
        x_pred1(:,t+1) = F*x_filt(:,t);
        V_pred1(:,:,t+1) = F*V_filt(:,:,t)*F'+Q;
    end
    Rhat = 0;
    for t=1:T
        Rhat = Rhat*(t-1)/t+(y(t)-H*x_pred1(:,t))^2/(H*V_pred1(:,:,t)*H'+R)/t;
    end
    ll = -T*log(Rhat)/2-T/2;
    for t=1:T
        ll = ll-log(H*V_pred1(:,:,t)*H'+R)/2;
    end
    mll = -ll;
end

