function [mll,g] = expGLM_ll(X,y,beta)
    mll = 0;
    g = zeros(size(X,2),1);
    for i=1:length(y)
        mll = mll+log(X(i,:)*beta)+y(i)/(X(i,:)*beta);
        g = g+X(i,:)'/(X(i,:)*beta)-y(i)*X(i,:)'/(X(i,:)*beta)^2;
    end
end

