function [FEVD, FEC] = DSGEforecasterror(theta,hMax)
%Function to construct forecast error covariance matrix and from one shock
%at a time at distance up to hMax

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);



%%%%%%%%%%%%%% Forecast Error Covariance Matrix
FEC = zeros(5,5);
for s = 1:hMax-1
    FEC = FEC + Phi1^s*Phi_eps*Phi_eps'*(Phi1^s)';
end

FEC = Psi1*FEC*Psi1';


%%%%%% Forecast error by individual shock
FEVD = zeros(4,4);
for j = 1:4
    I = zeros(4,4);
    I(j,j) = 1;
    for i = 1:4
        num = zeros(5,5);
        den = zeros(5,5);
        for s = 0:hMax-1
            num = num + Phi1^s*Phi_eps*I*Phi_eps'*(Phi1^s)';
            den = den + Phi1^s*(Phi_eps*Phi_eps')*(Phi1^s)';
        end
        num = Psi1*num*Psi1';
        den = Psi1*den*Psi1';
        FEVD(i,j) = num(i,i)./den(i,i);
    end
end

%FEVD(i,j) is contribution from j to i,

end