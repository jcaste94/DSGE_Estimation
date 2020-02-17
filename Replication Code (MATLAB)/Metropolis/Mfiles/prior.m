function [prior] = prior(Theta)

% The priors considered are:
%  1.) 100*(1/Beta-1) is GAMMA with mean 0.5 and std 0.5
%  2.) 100*log(Gamma) is NORMAL with mean .75 and std 0.5
%  3.) Lambda is GAMMA with mean .2 and std 0.2
%  4.) 100*log(PiStar) is GAMMA with mean 1 and std 0.5
%  5.) zetaP is BETA with mean 0.7 and std 0.15
%  6.) 1/(1+nu) is GAMMA with mean 1.50 and std 0.75
%  7.) rhoPhi is UNIFORM with mean 0 and std 1
%  8.) rhoLambda is UNIFORM with mean 0 and std 1
%  9.) rhoZ is UNIFORM with mean 0 and std 1
%  10.) 100*sigmaPhi is InvGamma with mean 2 and std 4.0
%  11.) 100*sigmaLambda is InvGamma with mean 0.5 and std 4.0
%  12.) 100*sigmaZ is InvGamma with mean 2 and std 4.0
%  13.) 100*sigmaR is InvGamma with mean 0.5 and std 4

% % % % NOTE: mean refers to para1, and std refers to para2

P = ones(13,1);

% % GAMMA pdf
para1 = [0.5,  0.2, 1, 1.50];
para2 = [0.5,  0.2,  0.5, 0.75];

b = para2.^2./para1;
a = para1./b;


X = [Theta(1), Theta(3), Theta(4), Theta(6)];

P([1,3,4, 6]) = gampdf(X,a,b);


% % NORMAL pdf 
para1 = 0.75;
para2 = 0.5;

P(2) = normpdf(Theta(2), para1, para2);

% % BETA pdf
para1 = 0.7;
para2 = 0.15;

a = (1-para1).*para1.^2./para2.^2 - para1;
b = a.*(1./para1 - 1);

P(5) = betapdf(Theta(5), a, b);

% % UNIFORM pdf
X = [Theta(7), Theta(8), Theta(9)];
if max( X<=0 ) ==1  || min(X<=1)==0
    P(7:9) = 0;
else
    P(7:9) = 1;
end

%use beta for now (when maximizing over prior)
% X = [Theta(7), Theta(8), Theta(9)];
% 
% para1 = .8*ones(1,3);
% para2 = .1*ones(1,3);
% 
% a = (1-para1).*para1.^2./para2.^2 - para1;
% b = a.*(1./para1 - 1);
% 
% P(7:9) = betapdf(Theta(7:9), a, b);



%Inverse gamma
para1 = [2, 0.5, 2, 0.5]';
para2 = [4, 4, 4, 4]';

X = Theta(10:13)';

pdf = exp(lnpdfig(X, para1, para2));

if max( X<=0)==1
    P(10:13) = 0;
else
    P(10:13) = pdf; 
end


f = prod(P);

prior =log(f);


% % Prior from inverse Gamma
% P8 = exp(lnpdfig(Theta(8),.01,4));
% P9 = exp(lnpdfig(Theta(9),.01,4));
% 
% f = P1*P2*P3*P4*P5*P6*P7*P8*P9;

% if f <= 10^(-13)
%    prior = -Inf;
% else
%    prior =log(f);
% end

%Theta

if prior==-Inf
    prior = -1e10;
end


end
function y = lnpdfig(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
%y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);

y = log(2) - gammaln(b./2) + (b./2).*log(b.*a.^2/2) - ( (b+1)/2 ).*log(x.^2) - b.*a.^2./(2*x.^2);

end




