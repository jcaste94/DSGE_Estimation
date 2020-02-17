% Readme for ANU lecture 
% - Estimation of DSGE Model
% 6/16/2013

%4/1/2015 - Jacob Warren
%NOTE: I changed the kalman function to reflect the linear model (deleted R and transposed A)
		     dsgeliki to use our data and model
		     prior to reflect our priors

There are 7 main files in the folder. 

1) Candidate.m                 : It computes the candidate density for Metropolis-Hastings algorith, 
                                 The output (MH_candidate.mat) is saved in the folder ``Matfiles".
	 
	%NOTE: RIGHT NOW, I GET ZERO GRADIENT FOR THE POSTERIOR MAXIMIZATION, SO I JUST USE THE ORIGINAL VALUES
	THE HESSIAN ALSO COMES BACK NON-POSITIVE DEFINITE, SO I USE IDENTITY

2) MetropolisHastings.m        : It runs the Metropolis-Hastings algorithm. It save the posterior draws 
                                 and logposterior density (mhdraws.mat) in the folder "Matfiles". The file 
                                 produces some diagnostic graph.

