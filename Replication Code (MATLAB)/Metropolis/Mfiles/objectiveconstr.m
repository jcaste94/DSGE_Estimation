function [objective] = objectiveconstr(Theta)
    
%    Theta(1) = exp(Theta(1));
%    Theta(2) = exp(Theta(2))/(1+exp(Theta(2)));
%    Theta(3) = exp(Theta(3));
%    Theta(4) = exp(Theta(4))/(1+exp(Theta(4)));
%    Theta(5) = exp(Theta(5))/(1+exp(Theta(5)));
%    Theta(8) = exp(Theta(8));
%    Theta(9) = exp(Theta(9));
    
prio = prior(Theta);
%Theta

save Matfiles/Theta

if prio==-Inf
    objective = -1e16;
else
    
liki = dsgeliki(Theta);

objective = (liki+prio);

end
objective = -objective;
% 
% disp('                                                                  ');    
% disp('                            THEta IS:                ');
% disp('                                                                  ');
% disp(num2str(Theta))
% disp('                                                                  ');            
% disp('                            objective IS:                ');
% disp('                                                                  ');
% disp(num2str(objective))
% disp('                                                                  ');   
%   
