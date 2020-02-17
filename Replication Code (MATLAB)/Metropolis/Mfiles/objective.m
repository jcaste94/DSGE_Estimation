function [objective] = objective(data,Theta,c)

global f fx fxp fy fyp
    
if c==1
    save parameter3 Theta
    Theta(1) = exp(Theta(1));
    Theta(2) = exp(Theta(2))/(1+exp(Theta(2)));
    Theta(3) = exp(Theta(3));
    Theta(4) = exp(Theta(4))/(1+exp(Theta(4)));
    Theta(5) = exp(Theta(5))/(1+exp(Theta(5)));
    Theta(8) = exp(Theta(8));
    Theta(9) = exp(Theta(9));
end

prio = prior(Theta);

if prio==-Inf
    objective = -1000000000000000;
else

liki = dsgeliki(data,Theta);

objective = liki+prio;
end

