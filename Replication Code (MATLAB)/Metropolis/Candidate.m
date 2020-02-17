%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%        Constructing the Candidate Density for MH Algorithm
%
%
%
% Author: Luigi Bocola         lbocola@sas.upenn.edu
% Date  : 06/16/2013
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

cd('/home/jascob/Dropbox/FileExchange_FS_JW/Matlab/Metropolis')


clear all
clc
close all
delete *.asv

tic

l = path;

path('Mfiles',path);
path('Optimization Routines',path);
path('LRE',path);
path('Matfiles',path);

disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: THE CANDIDATE DENSITY      ');
disp('                                                                  ');

%=========================================================================
%                          POSTERIOR MODE
%=========================================================================

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

%simulate data
[Yhat, s_true] = estimDSGE(theta, 500, 100);

save('Data/Yhat','Yhat')



%perform transformations:
%  1.) 100(1/Beta-1) 
%  2.) 100*log(Gamma)
%  3.) Lambda 
%  4.) 100*log(PiStar) 
%  5.) zetaP 
%  6.) 1/(1+nu) 
%  7.) rhoPhi 
%  8.) rhoLambda 
%  9.) rhoZ 
%  10.) 100*sigmaPhi 
%  11.) 100*sigmaLambda 
%  12.) 100*sigmaZ 
%  13.) 100*sigmaR 

theta(1) = 100*(1/theta(1) - 1);
theta(2) = 100*log(theta(2));
theta(4) = 100*log(theta(4));
theta(6) = 1/(1+theta(6));
theta(10:13) = 100*theta(10:13);


param = theta;


disp('                                                                  ');
disp('         *******STEP 1: RECOVERING THE POSTERIOR MODE....*********');
disp('                                                                  ');

objective = @(theta) prior(theta)*(-1) + dsgeliki(theta)*(-1);
[fh,x,gh,H,itct,fcount,retcodeh] = csminwel(objective,param,eye(length(param)),[] ,10^(-5),200);


Theta=x;

Theta(1) = 1/(1 + Theta(1)/100 );
Theta(2) = exp(Theta(2)/100);
Theta(4) = exp(Theta(4)/100);
Theta(6) = 1/Theta(6) - 1;
Theta(10:13) = Theta(10:13)/100;


disp('                                                                  ');    
disp('                            THE POSTERIOR MODE IS:                ');
disp('                                                                  ');
disp('   BETA      GAMMA      LAMBDA        Pi*        ZETA_p          NU    RHO_phi      RHO_lambda    RHO_Z  SIGMA_Phi  SIGMA_lambda  SIGMA_Z     SIGMA_R ');
disp(num2str(Theta))
disp('                                                                  ');                 
  

  
%=========================================================================
%                          CANDIDATE DENSITY
%=========================================================================

disp('                                                                  ');
disp('            *******STEP 2: HESSIAN AT SOLUTION....************'    );
disp('                                                                  ');


mode  = x; 

Sigma = nhess(objective,mode);
eig(Sigma)

Sigma = inv(Sigma);




%=========================================================================
%                            SAVE RESULTS
%=========================================================================

save Matfiles/MH_candidate Sigma mode 

path(l);

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;



