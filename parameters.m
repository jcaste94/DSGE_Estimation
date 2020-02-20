function theta = parameters

bbeta = 1/1.01;
ggamma = exp(.005);
llambda = .15;
piStar = exp(.005);
zetaP = .65;
nu = 0;
rhoPhi = .94;
rhoLambda = .88;
rhoZ = .13;
sigmaPhi = .01;
sigmaLambda = .01;
sigmaZ = .01;
sigmaR = .01;

%put into theta
theta(1) = bbeta;
theta(2) = ggamma;
theta(3) = llambda;
theta(4) = piStar;
theta(5) = zetaP;
theta(6) = nu;
theta(7) = rhoPhi;
theta(8) = rhoLambda;
theta(9) = rhoZ;
theta(10) = sigmaPhi;
theta(11) = sigmaLambda;
theta(12) = sigmaZ;
theta(13) = sigmaR;

end