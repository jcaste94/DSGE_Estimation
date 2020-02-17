function [Phi_star, Sigma_star] = DSGE_VAR_3(theta,p)
%Computed VAR approximation to DSGE model
%Note that default is to use the transformed parameters in
%DSGEautocovarIRF, but this can be changed if one wants to use a different
%set of tranformations to find the monetary policy shock.
%
%
%
% p is the order of VAR to calculate, theta is the set of DSGE model
% parameters


[~, ~, Psi0, ~] =  DSGE_soln_matrices(theta);   %get long-run mean
Psi0_temp = [Psi0(4) - Psi0(3)/theta(1); Psi0(1); Psi0(3)];
Psi0 = Psi0_temp;
m = length(Psi0);    %Number of variables in the VAR


gamma_0 = DSGEautocovarIRF(theta,0);

E_xx = zeros(m*p);
for i = 1:p
    E_xx((i-1)*m+1:i*m, (i-1)*m+1:i*m) = gamma_0 + Psi0*Psi0';
    
    for j = i:p
        h = abs(i-j);
        gamma_h = DSGEautocovarIRF(theta,h);
        
        E_xx( (i-1)*m+1:i*m , (j-1)*m + 1: j*m ) = gamma_h + Psi0*Psi0';
        E_xx( (j-1)*m + 1: j*m, (i-1)*m+1:i*m ) = gamma_h + Psi0*Psi0';
    end
end

E_yy = E_xx;
E_xx = [E_xx, repmat(Psi0,p,1); repmat(Psi0,p,1)', 1];


gamma_1 = DSGEautocovarIRF(theta,1);

E_xy = zeros(m*p);
for i = 1:p
    E_xy((i-1)*m+1:i*m, (i-1)*m+1:i*m) = gamma_1 + Psi0*Psi0';
    
    for j = i:p
        h = abs(i-j) + 1;
        gamma_h = DSGEautocovarIRF(theta,h);
        
        E_xy( (i-1)*m+1:i*m , (j-1)*m + 1: j*m ) = gamma_h + Psi0*Psi0';
        E_xy( (j-1)*m + 1: j*m, (i-1)*m+1:i*m ) = gamma_h + Psi0*Psi0';
    end
end


E_xy = [E_xy; repmat(Psi0,p,1)'];

Phi_star = inv(E_xx)*E_xy;
Phi_star(:, m+1:end) = [eye(m*(p-1)); zeros(m + 1, m*(p-1))];

Sigma_star = E_yy - E_xy'*inv(E_xx)*E_xy;
Sigma_star = Sigma_star(1:m,1:m);

end