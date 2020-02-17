function [s_ss, s_yy, f_ss, f_yy] = DSGEspectrum(theta, omega)
%compute the spectrum for the DSGE model at a vector of frequencies, omega

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);

%%%%%%%%%%%%% Compute spectrum at frequencies omega
i = sqrt(-1);
s_ss = zeros(5,5,length(omega));
for o = 1:length(omega)
    Sl = eye(5,5) - Phi1'*exp(i*omega(o));
    s_ss(:,:,o) = 1/(2*pi)*inv(Sl)*Phi_eps*Phi_eps'*inv(Sl');   %this is from the gauss code
end

%%%%%%%%%%%%% Spectral Density for observables
s_yy = zeros(4,4,length(omega));

for o = 1:length(omega)
    s_yy(:,:,o) = Psi1*s_ss(:,:,o)*Psi1';
end


%%%%%%%%%%%%% Contribution of shock j to spectral density:
i = sqrt(-1);
f_ss = zeros(5,5,length(omega));
f_yy = zeros(4,5,length(omega));
for j = 1:4
    I = zeros(4,4);
    I(j,j) = 1;
    for o = 1:length(omega)
        Sl = eye(5,5) - Phi1'*exp(i*omega(o));
        f_ss(:,:,o) = real( 1/(2*pi)*inv(Sl)*Phi_eps*I*Phi_eps'*inv(Sl') );   %this is from the gauss code
        
        %turn into observables
        f_yy(:,j,o) = diag(Psi1*f_ss(:,:,o)*Psi1');
    end
end


end