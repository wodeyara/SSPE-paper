function [x_new,P_new] = oneStepKFupdate_sspp(x,y,phi,M,Q,R,P)
% one step kalman filter update given the current time point's x, the
% observed y, the variance paramter estimates of Q, R and P (For x itself)

x_one = phi *x;
P_one = phi *P * phi' + Q;

K_one = (P_one*M) /(M'*P_one*M + R); 

x_new = x_one + K_one * (y - M'*x_one);
P_new = P_one - K_one * M' * P_one;
