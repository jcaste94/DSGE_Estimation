function [Gamma_ss_h, Gamma_yy_h, Gamma_ys_h] = DSGEautocovar(theta, h)
%Compute autocovariance function for lag length h at parameters theta


[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);

%%%%%% Solve Matrix Ricatti Equation to get autocovariance of s_t at 0
Phi_eps2 = Phi_eps*Phi_eps';
vecGamma_ss_0 = inv(eye(5*5,5*5) - kron(Phi1,Phi1))*Phi_eps2(:);
Gamma_ss_0 = reshape(vecGamma_ss_0,5,5);

%%%%% Construct autocovariance at desired lag
if(h > 0)
    Gamma_ss_h = Phi1^h*Gamma_ss_0;
elseif(h<0)
     Gamma_ss_h = (Phi1^abs(h)*Gamma_ss_0)';
else
    Gamma_ss_h = Gamma_ss_0;
end
Gamma_yy_h = Psi1*Gamma_ss_h*Psi1';
Gamma_ys_h = Psi1*Gamma_ss_h;

end