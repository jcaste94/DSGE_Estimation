function irf = ImpulseResponse(theta,hMax, from)
%compute impulse response functions for s

%from is an integer indicating which variable originates the shock

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);


%%%%%%% Compute impulse response function

if hMax==0
    irf = Psi1*Phi_eps(:,from);
else
    irf = Psi1*(Phi1^hMax)*Phi_eps(:,from);
end

end
