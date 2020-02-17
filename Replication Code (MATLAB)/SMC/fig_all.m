% ----------------------------------------------------------------
% Figure file for SMC example
% - This file creates/saves figures based on outputs from main_ss.m
% Updates
% 5/13/2014: scatter plot 
% ----------------------------------------------------------------
% Minchul Shin [mincshin@sas.upenn.edu]
% Last modified: 5/13/2014
% ----------------------------------------------------------------
clc; clear all; close all;

% add the toolbox
addpath([pwd, '\toolbox_ss']);

%% General Setting - Specificaitons/Path
modelname = 'HS'; % Herbst and Schorfheide
priorname = '0';
specname  = [modelname, '_p', priorname]; % prior 
smcrun    = '1'; % chain number
propmode  = 1;   % 1-RWMH, 2-Mixture MH
T0_       = 1;   % starting point of likelihood evaluation

rstname = [specname, '_mode',num2str(propmode),'_run',smcrun]; 

% path
workpath = pwd;
priopath = workpath; %prior path

currentFolder = pwd;
lastSlashPosition = find(currentFolder == '/', 1, 'last');
parentFolder = currentFolder(1:lastSlashPosition-1);

savepath = [parentFolder, '/figuresv2/smc/'];
chk_dir(savepath);

%% Load results
loadfilename = [specname, '_mode',num2str(propmode),'_run',smcrun,'.mat'];
cd(savepath);
load(loadfilename);
cd(workpath);

%% Figure 1: Waterfall plot


 pnames = strvcat('\beta','\gamma','\lambda','\pi^{*}','\zeta_{p}', '\nu',...
     '\rho_{\phi}','\rho_{\lambda}','\rho_{Z}','\sigma_{\phi}', ...
     '\sigma_{\lambda}','\sigma_{Z}','\sigma_R');
 
 
 for parai = [5, 10]
     fig = figure(1);
     setmyfig(fig, [1.7, 1.2, 8, 7], version)
     
     %parai = 1; %parameter index
     
     itsel = 1:1:tune.nphi;
     nsel  = length(itsel);
     
     if parai == 5
       bins = 0:0.01:1.2;
     else
         bins = 0.5:.01:4;
     end
     post = zeros(nsel, length(bins));
     for i = 1:nsel
         
         phisel  = itsel(i);
         para    = squeeze(parasim(phisel, :, parai));
         [id] = multinomial_resampling(wtsim(:, phisel)');
         para    = para(id);
         
         post(i, :) = ksdensity(para, bins);
     end
     
     colscale = (0.65:-0.05:0);
     colormap(repmat(colscale', 1, 3));
          
     phigrid = (1:1:tune.nphi)';
     waterfall(bins, phigrid, post)
     
     xlim([min(bins), max(bins)])
     ylim([0, tune.nphi])

     set(gca,'fontsize',20,'fontweight','demi')
     xlabel(['$',pnames(parai,:),'$'], 'fontsize', 20,'interpreter', 'latex')
     ylabel(['$N_{\phi}$'], 'fontsize', 20, 'interpreter', 'latex')
     
     cd(savepath);
     savefilename = ['fig_',rstname, '_waterfall_', num2str(parai), '.png'];
     saveas(fig, savefilename);
     cd(workpath);
     
 end


