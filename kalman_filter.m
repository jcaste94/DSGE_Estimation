function [L,hi,lo,postMean,postVar] = kalman_filter(Y,Psi0,...
    Psi1, Sigma_u, Phi1,Phi_eps, s0Hat,alpha, P)
%function outputs the likelihood value (L), the upper and lower bounds (hi
%and lo) of an "alpha" percentage credible set, and E(s_t+1 | Y_t) (s), in
%addition to the posterior mean and variance of the state|Y


%forecasting step
sHat_t = Phi1*s0Hat;
P_t = Phi1'*P*Phi1 + Phi_eps*Phi_eps';

%predicting Y:
yHat_t = Psi0 + (Psi1*sHat_t)';
F_t = Psi1*P_t*Psi1' + Sigma_u;
v = Y-yHat_t; 

%compute likelihood
L = mvnpdf(v,0,F_t);

%update using Bayes Theorem (updated values are the prior for next period's
%iteration of the filter as well)
postMean = sHat_t + P_t*Psi1'*inv(F_t)*v';
postVar = P_t - P_t*Psi1'*inv(F_t)*Psi1*P_t;

%impose zero variance for lagged variable
postVar(5,:) = zeros(5,1);
postVar(:,5) = zeros(5,1);

%estimate alpha-level confidence interval
lo = postMean(1:4) - norminv(1 - alpha/2, 0, sqrt(diag(postVar(1:4,1:4))));
hi = postMean(1:4) + norminv(1 - alpha/2, 0, sqrt(diag(postVar(1:4,1:4))));  

end


