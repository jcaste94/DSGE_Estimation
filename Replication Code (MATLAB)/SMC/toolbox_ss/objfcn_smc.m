function [lnpost, lnpY, lnprio, error] = objfcn_smc(para, phi_smc)
% objfcn_smc
%  - phi_smc is tempering parameter
%  - evaluates "TEMPERED" posterior density given parameter bounds
%  - Modified based on professor Schorfheide's gauss code
%  - by Minchul Shin (UPenn)
% updated: 2014-05-08


%check prior to make sure it is in bounds
prior_val = prior(para);


if prior_val == -Inf
    lnpost  = -inf;
    lnpY    = -inf;
    lnprio = -inf;
    error   = 10;
else
     % evaluation of likelihood
     
     lnpY = dsgeliki(para);

    % evaluate the prior distribution
    lnprio = prior_val;
    
    % log posterior
    lnpost = (phi_smc*lnpY + lnprio);
end







