%==========================================================================
%                       DSGE MODEL ESTIMATION
%                   Metropolis-Hastings Algorithm
%
%
% Author: Luigi Bocola         lbocola@sas.upenn.edu
% Date  : 06/16/2013
%
% Edited: Jacob Warren         jacobwar@sas.upenn.edu
% Date  : 11/16/2015
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================
clc
clear all
close all
delete *.asv

tic

l = path;

path('Mfiles',path);
path('Optimization Routines',path);
path('LRE',path);
path('Matfiles',path);

currentFolder = pwd;
lastSlashPosition = find(currentFolder == '/', 1, 'last');
parentFolder = currentFolder(1:lastSlashPosition-1);

savePath = [parentFolder, '/figuresv2/metropolis/'];

disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: METROPOLIS-HASTINGS        ');
disp('                                                                  ');

%=========================================================================
%                  METROPOLIS-HASTINGS ALGORITHM
% (Report the Acceptance Rate and Recursive Averages Every 1000 draws)
%=========================================================================

load MH_candidate

Nsim          = 30000;
disp('                                                                  ');
disp('                                                                  ');

rng(10)

c             = 0.075;
Nburn         = int32(0.25*Nsim);
Nsim = Nsim+Nburn;
Thetasim      = zeros(Nsim,length(mode));
Thetasim(1,:) = draw_prior;             %make initial draw from prior
%Thetasim(1,:) = mode;            %make initial draw from mode
accept        = 0;
obj           = dsgeliki(Thetasim(1,:)) + prior(Thetasim(1,:));
counter       = 0;
logposterior  = obj*ones(Nsim,1);

for i=1:Nsim
    
    Thetac = mvnrnd(Thetasim(i,:),c*Sigma);
    
    prioc=prior(Thetac);
    
    if prioc== -1e10
        Thetasim(i+1,:) = Thetasim(i,:);
    else
                
        likic = dsgeliki(Thetac);
        objc  = prioc+likic;

        alpha = min(1,exp(objc-obj));
        
        u = rand(1);
        
        if u<=alpha
            Thetasim(i+1,:)   = Thetac;
            accept            = accept+1;
            obj               = objc;
            logposterior(i+1) = objc;
        else
            
            Thetasim(i+1,:)   = Thetasim(i,:);
            logposterior(i+1) = obj;
            
        end
        
        acceptancerate     = accept/double(i);
        counter            = counter + 1;
        
        
        if counter==500
            disp('                                                                  ');
            disp(['                               DRAW NUMBER:', num2str(i)]         );
            disp('                                                                  ');
            disp('                                                                  ');
            disp(['                           ACCEPTANCE RATE:', num2str(acceptancerate)]);
            disp('                                                                  ');
            disp('                                                                  ');
            disp('                            RECURSIVE AVERAGES                    ');
            disp('                                                                  ');
            disp(num2str(mean(Thetasim(1:i,:))));
            disp('                                                                  ');
            
            counter = 0;
        end
    end
end

%acceptance rate is 

%%%%%%%%%%%%%%%%%% Print the raw draws from the MCMC %%%%%%%%%%%%%%%%%%%%
[Nsim,Npam] = size(Thetasim);

 pnames = strvcat('\beta','\gamma','\lambda','\pi^{*}','\zeta_{p}', '\nu',...
     '\rho_{\phi}','\rho_{\lambda}','\rho_{Z}','\sigma_{\phi}', ...
     '\sigma_{\lambda}','\sigma_{Z}','\sigma_R');
 


figure('Position',[20,20,900,600],'Name',...
    'MCMC draws','Color','w')

for i=1:Npam
    
    subplot(ceil(Npam/3),3,i), plot(Thetasim(:,i),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
    title(pnames(i,:),'FontSize',12,'FontWeight','bold');
end

close all

j = 5;
plot(Thetasim(:,j),'LineStyle','-','Color','b','LineWidth',4)
set(gca,'fontsize',20,'fontweight','demi')
box off
print('-dpng', [savePath 'draw_zetaP'])

j = 10;
plot(Thetasim(:,j),'LineStyle','-','Color','b','LineWidth',4)
set(gca,'fontsize',20,'fontweight','demi')
box off
print('-dpng', [savePath 'draw_sigma_phi'])





Thetasim    = Thetasim(Nburn:end,:);

logposterior= logposterior(Nburn:end);

save Matfiles/mhdraws Thetasim logposterior   % Save posterior draws

load Matfiles/mhdraws Thetasim logposterior
[Nsim,Npam] = size(Thetasim);

%=========================================================================
%                  FIGURE 1: RECURSIVE AVERAGES
%=========================================================================



figure('Position',[20,20,900,600],'Name',...
    'Recursive Averages','Color','w')

rmean = zeros(Nsim,Npam);


for j=1:Nsim
    rmean(j,:) = mean(Thetasim(1:j,:),1);
end

for j=1:Npam
    
    subplot(ceil(Npam/3),3,j), plot(rmean(:,j),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
    title(pnames(j,:),'FontSize',12,'FontWeight','bold');
end
hold off

close all
%plot zetaP and sigma_phi
j = 5;
plot(Nburn:length(rmean)+Nburn-1,rmean(:,j),'LineStyle','-','Color','b','LineWidth',4)
set(gca,'fontsize',20,'fontweight','demi')
box off
ylim([.6, .7])
print('-dpng', [savePath 'recursive_mean_zetaP'])

j = 10;
plot(Nburn:length(rmean)+Nburn-1,rmean(:,j),'LineStyle','-','Color','b','LineWidth',4)
set(gca,'fontsize',20,'fontweight','demi')
box off
ylim([.9, .95])
print('-dpng', [savePath 'recursive_mean_sigma_phi'])





%=========================================================================
%                  FIGURE 2: POSTERIOR MARGINAL DENSITIES
%=========================================================================


figure('Position',[20,20,900,600],'Name',...
    'Posterior Marginal Densities','Color','w')


for i=1:Npam
    xmin = min(Thetasim(:,i));
    xmax = max(Thetasim(:,i));
    grid = linspace(xmin,xmax,100);
    u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
    subplot(ceil(Npam/3),3,i), plot(grid,ksdensity(Thetasim(:,i)),'LineStyle','-','Color','b',...
        'LineWidth',2.5), hold on
    plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
        'Color','black','LineWidth',2.5 ), hold off
    axis([xmin xmax 0 u]);
    title(pnames(i,:),'FontSize',12,'FontWeight','bold');
end
hold off


numDraws = 5000;
priorsim = zeros(numDraws, size(Thetasim,2));
for j = 1:numDraws
    priorsim(j,:) = draw_prior;
end

close all

i = 5;
xmin = .5;
xmax = 1;
grid = linspace(xmin,xmax,100);
u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
plot(grid,ksdensity(Thetasim(:,i), grid),'LineStyle','-','Color','b',...
    'LineWidth',4), hold on
plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
    'Color',.5 + [0,0,0],'LineWidth',2.5 ), 
plot(grid, ksdensity(priorsim(:,i), grid), 'LineStyle', '--', 'Color','black',...
    'LineWidth',4)
hold off
axis([xmin xmax 0 u]);
box off
set(gca,'fontsize',20,'fontweight','demi')
print('-dpng', [savePath 'density_zetaP'])

i = 10;
xmin = .4;
xmax = 1.2;
grid = linspace(xmin,xmax,100);
u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
plot(grid,ksdensity(Thetasim(:,i), grid),'LineStyle','-','Color','b',...
    'LineWidth',4), hold on
plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
    'Color',.5 + [0,0,0],'LineWidth',2.5 ), 
plot(grid, ksdensity(priorsim(:,i),grid), 'LineStyle', '--', 'Color','black',...
    'LineWidth',4)
hold off
axis([xmin xmax 0 u]);
box off
set(gca,'fontsize',20,'fontweight','demi')

print('-dpng', [savePath 'density_sigma_phi'])


disp('                                                                  ');
disp(['                     ELAPSED TIME:   ', num2str(toc)]             );

elapsedtime=toc;

path(l);


