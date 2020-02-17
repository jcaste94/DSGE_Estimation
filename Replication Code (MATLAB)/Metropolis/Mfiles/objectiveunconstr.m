function [objective] = objectiveunconstr(Theta)
    
prio = prior(Theta);


if prio==-Inf
    objective = -1000000000000000;
else

liki = dsgeliki(Theta);

%load Data/FedFunds
%load Data/GDP
%load Data/GDPdef
%load Data/Labor

%Yhat = [GDPdat LabShare GDPDEF FedFundsQuarterly];
 

%Sigma_u = zeros(4);
%alpha = .1;
%[s_mean, P_mean, s_hi, s_lo, L] = kal_wrapper(Theta,Yhat,alpha, Sigma_u);

%liki = sum(log(L));


objective = -(liki+prio);

end

