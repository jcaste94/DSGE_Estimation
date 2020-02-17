%--------------------------------------------
%--------------------------------------------
% File to perform MLE of DSGE model
%--------------------------------------------
% Jacob Warren
% jacobwar@sas.upenn.edu
%--------------------------------------------
% 11/19/2015
%--------------------------------------------
%--------------------------------------------

clc, clear all, close all

%add helper functions:
addpath([pwd, '/helper_functions']);

%path to save pictures
picPath = [pwd, '/figuresv2/MLE/'];

% pass in parameters from paper
theta = parameters;


%% Sample size 80

rng(101)                                   % Set the seed for reproducibility
numReps = 100;                             % Number of simulations
sampleSizes = [80,200];                    %Sample sizes to consider
zeta_opt = zeros(numReps,length(sampleSizes));    %matrix to store optimal zetas
Tau = linspace(.45,.95, 50);   %% vary zetaP across the grid


parfor n = 1:length(sampleSizes)
    N = sampleSizes(n);
    
    rng(n*100)
    
    for jj = 1:numReps
        
        [Yhat, s_true] = estimDSGE(theta, N, 100);
        
       
        Q = zeros(length(Tau),1);    %objective function
        theta_new = theta;
        for z = 1:length(Tau)    %discrete grid for evaluating the likelihood
            zeta = Tau(z);
            theta_new(5) = zeta;
            
            %kalman filter for likelihood
            Sigma_u = zeros(4);
            alpha = .1;
            [~,~,~,~, liki] = kal_wrapper(theta_new, Yhat, alpha, Sigma_u);

            
            Q(z) = sum(log(liki));
        end
        
        [~,minLoc] = max(Q);
        
        zeta_opt(jj,n) = Tau(minLoc);
        
    end
end
zetaP = theta(5);

[d, grid] = ksdensity(zeta_opt(:,1));
plot(grid,d, 'linewidth',4, 'Color','blue', 'linestyle',':')
hold on

[d, grid] = ksdensity(zeta_opt(:,2));
plot(grid,d, 'linewidth',4, 'Color','black', 'linestyle','--')

plot([zetaP zetaP], [0 max(d)+10],...
    'Color','red','LineWidth',2.5 ), hold off
ylim([0,max(d)+1])
xlim([0,1])
set(gca, 'Xtick',0:.2:1)
set(gca,'fontsize',20,'fontweight','demi')
box off
print('-dpng',[picPath, 'samp_distr_mle'])

%Blue (dotted) is sample size 80, black (dashed) is sample size 200, red
%(solid) is true value


%plot single likelihood:

[Yhat, s_true] = estimDSGE(theta, 200, 100);

Q = zeros(length(Tau),1);    %objective function
theta_new = theta;
Tau2 = linspace(.45,.75, length(Tau));
for z = 1:length(Tau2)    %discrete grid for evaluating the likelihood
    zeta = Tau2(z);
    theta_new(5) = zeta;
    
    %kalman filter for likelihood
    Sigma_u = zeros(4);

    alpha = .1;
    [~,~,~,~, liki] = kal_wrapper(theta_new, Yhat, alpha, Sigma_u);
    
    Q(z) = sum(log(liki));
end
[a,b] = max(Q);
Tau2(b)


plot([zetaP zetaP], [0 max(Q)+20],...
    'Color','red','LineWidth',2.5 )
hold on
plot(Tau2, Q, 'linewidth',4, 'Color','blue')

ylim([min(Q)-20,max(Q)+20])
set(gca, 'Xtick',0:.1:1)
set(gca,'fontsize',20,'fontweight','demi')

box off
print('-dpng',[picPath, 'likelihood'])

