function [liki] = dsgeliki(para)

% This function computes the likelihood of the DSGE model.
% Input = para: Vector of Structural Parameters
%         constr: Indicator equal to 1 if parameters are constrained
% Output= Likelihood Function

%para(3) = exp((1/4)*para(3));

%[T1, TC, T0, TETA, RC, retcode] = model_solution(para);

%if retcode==0
%    
%    load data
%
%    data(:,1) = data(:,1)-log(para(3))*linspace(1,size(data,1),size(data,1))';
%    
%    [A,B,H,R,Se,Phi] = sysmat(T1,T0,para);
%    
%    liki = kalman(A,B,H,R,Se,Phi,data);
%    
%    liki = sum(liki);
%    
%else
%    
%    liki=-1000000000000;
%    
%end


load Data/Yhat

       

%reverse transformations in para:
para(1) = 1/(1 + para(1)/100 );
para(2) = exp(para(2)/100);
para(4) = exp(para(4)/100);
para(6) = 1/para(6) - 1;
para(10:13) = para(10:13)/100;


 
[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(para);


%if retcode==0
    
%Sigma_u = ones(4)*1e-5;
Sigma_u = zeros(4);
 %  liki = kalman(Psi0 ,Psi1, zeros(4), eye(5,5) , Phi_eps*Phi_eps',Phi1,Yhat);
   [s_mean, P_mean, s_hi, s_lo, Likelihood] = kal_wrapper(para,Yhat,.1, Sigma_u);

 %  liki1 = sum(liki(75:end));
 %  liki2 = sum(Likelihood(75:end));
   
    
 %   liki = sum(liki);
%    liki = sum(real(liki));
      liki = sum(Likelihood);
%      liki
% 
%       if isreal(liki) == 1
%           liki = sum(liki);
%       else
%           liki = sum(Likelihood);
%       end
%     
%else
     
  %  liki=-1000000000000;
    
%end





end



