function [s_mean, P_mean, s_hi, s_lo, Likelihood] = kal_wrapper(theta,Y,alpha, Sigma_u)

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);

rhoPhi = theta(7);
rhoLambda = theta(8);
rhoZ = theta(9);
sigmaPhi = theta(10);
sigmaLambda = theta(11);
sigmaZ = theta(12);
sigmaR = theta(13);

%find out the number of series that are given, and truncate Psi matrices if
%necessary

p = size(Y,2);
Psi0 = Psi0(1:p);
Psi1 = Psi1(1:p,1:5);
Sigma_u = Sigma_u(1:p,1:p);

%extract the coefficients for x_{t-1}, to construct confidence intervals
Phi_past = Phi1(5,1:4);

%initial values

P0 = diag([sigmaPhi^2/(1-rhoPhi^2), sigmaLambda^2/(1-rhoLambda^2), sigmaZ^2/(1-rhoZ^2), sigmaR^2, 0]);

s0 = zeros(5,1);
P = P0;

%initialize objects to store returns in
Likelihood = zeros(length(Y),1);
s_hi = zeros(length(Y),4);
s_lo = zeros(length(Y),4);
s_mean = zeros(length(Y),5);
P_mean = zeros(length(Y),5,5);

%perform the filtering
for t = 1:length(Y)
    %perform filtering at date t
    [L,hi,lo,postMean,postVar] =  kalman_filter(Y(t,:),Psi0, Psi1, Sigma_u, ...
        Phi1,Phi_eps, s0 ,alpha, P);
    
    %store output values
    Likelihood(t) = L;
    s_hi(t,1:4) = hi;
    s_lo(t,1:4) = lo;
    s_mean(t,:) = postMean;
    P_mean(t,:,:) = postVar;
    P = postVar;
    
    %compute confidence intervals for x_{t-1}
    if t>1
        s_lo(t,5) =  Phi_past*s_lo(t-1,1:4)';
        s_hi(t,5) = Phi_past*s_hi(t-1,1:4)';
    end
    
    s0 = postMean;
end

end