%--------------------------------------------
%--------------------------------------------
% File to perform moment matching on DSGE model
%--------------------------------------------
% Jacob Warren
% jacobwar@sas.upenn.edu
%--------------------------------------------
% 5/22/2015
%--------------------------------------------
%--------------------------------------------

clc, clear all, close all

%add helper functions:
addpath([pwd, '/helper_functions']);

%path to save pictures
picPath = [pwd, '/figuresv2/methodofmoments/'];

% pass in parameters from paper
theta = parameters;


%% Moment Matching

rng(101)                                   % Set the seed for reproducibility
numReps = 100;                             % Number of simulations
sampleSizes = [80,200];                     %Sample sizes to consider
zeta_opt = zeros(numReps,length(sampleSizes));    %matrix to store optimal zetas
zeta_opt_dsge = zeros(numReps,length(sampleSizes));    %matrix to store optimal zetas
Tau = linspace(.2,.95, 50);   %% vary zetaP across the grid
N_E = 100;          % Number of datasets to compute expected moments


for n = 1:length(sampleSizes)
    N = sampleSizes(n);
    
    
    parfor jj = 1:numReps
        rng(jj*N)

        [Yhat, s_true] = estimDSGE(theta, N, 100);        
        Y_true = Yhat(:,[1,3]);               %match only on output and inflation

        p = 2;      %VAR(2)
        
        Yp = Y_true(1:length(Y_true)-p,:);
        Xp = Y_true(2:length(Y_true)-(p-1),:);
        for j = 2:p
            Xp = [Xp, Y_true(j+1:length(Y_true)-(p-j),:)];
        end
        Xp = [ones(length(Xp),1),Xp];
        
        PhiHat_true = (Xp'*Xp)\(Xp'*Yp);
        sigmaHat = 1/length(Yp)*(Yp - Xp*PhiHat_true)'*(Yp - Xp*PhiHat_true);
        
         W = kron(inv(sigmaHat),Xp'*Xp);    %use optimal weight matrix      
       
        Q = zeros(length(Tau),1);
        
        theta_new = theta;
        for z = 1:length(Tau)
            zeta = Tau(z);
            theta_new(5) = zeta;
            
            %simulate new data
            mHat = zeros(size(PhiHat_true));
            
            %use the same random numbers for each simulation to decrease
            %monte carlo error
            rng(jj + 1e5)
            
            for n_e = 1:N_E
                [Yhat, s_true] = estimDSGE(theta_new, N, 100);
                
                Y = Yhat(:,[1,3]);
                
                Yp = Y(1:length(Y)-p,:);
                Xp = Y(2:length(Y)-(p-1),:);
                for j = 2:p
                    Xp = [Xp, Y(j+1:length(Y)-(p-j),:)];
                end
                Xp = [ones(length(Xp),1),Xp];
                
                PhiHat = (Xp'*Xp)\(Xp'*Yp);
                
                mHat = mHat + 1/N_E * PhiHat;
            end
            
            err = PhiHat_true - mHat;
            Q(z) = err(:)' * W * err(:);
        end
        
        [~,minLoc] = min(Q);
        
        zeta_opt(jj,n) = Tau(minLoc);
        
        % DSGE implied VAR coefficient
        theta_new = theta;
        for z = 1:length(Tau)
            zeta = Tau(z);
            theta_new(5) = zeta;
            
            [phiHat_dsge, ~, ~] = DSGE_VAR_2(theta_new, 2);
            
            mHat = phiHat_dsge(:,1:2);

            err = PhiHat_true - mHat;
            Q(z) = err(:)' * W * err(:);
        end
        
        [~,minLoc] = min(Q);
        
        zeta_opt_dsge(jj,n) = Tau(minLoc);
        
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
ylim([0,8])
xlim([0,1])
set(gca, 'Xtick',0:.2:1)
set(gca,'fontsize',20,'fontweight','demi')
box off
print('-dpng',[picPath, 'MM'])

%Blue (dotted) is sample size 80, black (dashed) is sample size 200, red
%(solid) is true value


hold off
figure(2)
[d, grid] = ksdensity(zeta_opt_dsge(:,1));
plot(grid,d, 'linewidth',4, 'Color','blue', 'linestyle',':')
hold on

[d, grid] = ksdensity(zeta_opt_dsge(:,2));
plot(grid,d, 'linewidth',4, 'Color','black', 'linestyle','--')

plot([zetaP zetaP], [0 max(d)+10],...
        'Color','red','LineWidth',2.5 ), hold off
ylim([0,8])
xlim([0,1])
set(gca,'fontsize',20,'fontweight','demi')
set(gca, 'Xtick',0:.2:1)
box off
print('-dpng',[picPath, 'MM_dsge'])

%Blue (dotted) is sample size 80, black (dashed) is sample size 200, red
%(solid) is true value
