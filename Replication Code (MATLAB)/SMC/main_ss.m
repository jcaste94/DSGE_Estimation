% ----------------------------------------------------------------
% Main file for SMC example
% ----------------------------------------------------------------
% - This file replicates Figure 1 in Herbst and Schorfheide (2014).
% - Code is based on Ed Herbst's code [http://edherbst.net/]. 
% ***If you use this code, please cite Herbst and Schorfheide (2014).
% ----------------------------------------------------------------
% Updates
% 5/13/2014: We don't use dlyap.m
% ----------------------------------------------------------------
% Minchul Shin [mincshin@sas.upenn.edu]
% Last modified: 5/13/2014
% ----------------------------------------------------------------
clc; clear all; close all;

cd('/home/jascob/Dropbox/FileExchange_FS_JW/Matlab/SMC')

%turn of parfor warnings
warning off

% add the toolbox
addpath([pwd, '/toolbox_ss']);

%% General Setting - Specificaitons/Path

propmode  = 1;   % 1-RWMH, 2-Mixture MH                   

% path
workpath = pwd;
%priopath = workpath; %prior path
savepath = [workpath, '/results_',specname];
chk_dir(savepath);

%% load data

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
[Yhat, s_true] = estimDSGE(theta, 100, 100);

save('Data/Yhat','Yhat')

YY = Yhat;

%Perform transformations

theta(1) = 100*(1/theta(1) - 1);
theta(2) = 100*log(theta(2));
theta(4) = 100*log(theta(4));
theta(6) = 1/(1+theta(6));
theta(10:13) = 100*theta(10:13);



%% SMC Setting - setting up a tuning variable, "tune"

% General
tune.npara = length(theta);      % # of parameters
tune.npart = 2*1024; % # of particles
tune.nphi  = 500;     % # of stage

%tune.npart = 300; % # of particles
%tune.nphi  = 10;     % # of stage

tune.lam   = 3;      % # bending coeff, lam = 1 means linear cooling schedule

% Create the tempering schedule
tune.phi = (1:1:tune.nphi)';
tune.phi = ((tune.phi-1)/(tune.nphi-1)).^tune.lam;

% tuning for MH algorithms
tune.c    = 0.1;                  % initial scale cov
tune.R    = 0.01*eye(tune.npara); % initial cov
tune.acpt = 0.25;                 % initial acpt rate
tune.trgt = 0.25;                 % target acpt rate
tune.alp  = 0.9;  % Mixture weight for mixture proposal

% run script
ss_smc; %this file will save parasim, wtsim, tune, others


