
%--------------------------------------------
%--------------------------------------------
% File to perform IRF matching on DSGE model
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
picPath = [pwd, '/figuresv2/irfmatching/'];

% pass in parameters from paper
theta = parameters;


%% Figure 1 
% plot DSGE IRF with the VAR approximation (population moments derived from
% autocovariance function)
% Compare DSGE IRF with VAR(1) 

%dsge implied VAR impulse response
hMax = 6;
p = 1;   %VAR(1)
[Phi_star, Sigma_star] = DSGE_VAR_3(theta,p);   %VAR(1)

Sigma_chol = chol(Sigma_star)';
F = Phi_star(1:length(Phi_star)-1,:)';   %companion Form

m = size(Sigma_chol,1);
M_n =[eye(m), zeros(m,(p-1)*m)];

C_hat = @(h) M_n*F^h*M_n';

VAR_IRF = zeros(hMax+1,m);
for h = 0:hMax
    VAR_IRF(h+1,:) = C_hat(h) * Sigma_chol(:,1);
end
%only use output and inflation (2 is output and 3 is inflation)


% Actual IRF from the dsge model
dsge_irf = zeros(4,hMax+1,5);
for j = 1:4
    for h = 0:hMax
        dsge_irf(:,h+1,j) = ImpulseResponse(theta,h,j);
    end
end

plot(0:hMax,cumsum(VAR_IRF(:,2))*100, 'color','blue', 'linewidth',4)
hold on
zzz = squeeze(dsge_irf(1,:,4));
plot(0:hMax,cumsum(zzz*100), 'linewidth',4, 'color','black', 'linestyle','--')
ylim([-.9, .1])
xlim([0,6])
xlabel('Horizon h','fontsize',22)
set(gca,'fontsize',20,'fontweight','demi')
box off
hold off
print('-dpng',[picPath, 'figure1_1'])

%dashed (black) line is DSGE IRF, solid (blue) is VAR implied


plot(0:hMax,VAR_IRF(:,3)*100, 'color','blue', 'linewidth',4)
hold on
zzz = squeeze(dsge_irf(3,:,4));
plot(0:hMax,zzz*100, 'linewidth',4, 'color','black','linestyle','--')
ylim([-.9, .1])
xlim([0,6])
xlabel('Horizon h','fontsize',22)
set(gca,'fontsize',20,'fontweight','demi')
box off
hold off
print('-dpng',[picPath, 'figure1_3'])

%% Sensitivity of IRF to zeta_p
zeta_other = 0.5;

theta_new  = theta;
theta_new(5) = zeta_other;

% Actual IRF from the dsge model
dsge_irf_new = zeros(4,hMax+1,5);
for j = 1:4
    for h = 0:hMax
        dsge_irf_new(:,h+1,j) = ImpulseResponse(theta_new,h,j);
    end
end


plot(0:hMax,cumsum(VAR_IRF(:,2))*100, 'color','blue', 'linewidth',4)
hold on
zzz = squeeze(dsge_irf(1,:,4));
plot(0:hMax,cumsum(zzz*100), 'linewidth',4, 'color','black', 'linestyle','--')
zzz = squeeze(dsge_irf_new(1,:,4));
plot(0:hMax,cumsum(zzz*100), 'linewidth',4, 'color','red', 'linestyle',':')
ylim([-.9, .1])
xlim([0,6])
xlabel('Horizon h','fontsize',22)
set(gca,'fontsize',20,'fontweight','demi')
box off
hold off

print('-dpng',[picPath, 'figure2_1'])

%dashed (black) line is DSGE IRF (zeta = .65), solid (blue) is VAR implied,
%dotted (red) lines is DSGE (zeta = .5)


plot(0:hMax,VAR_IRF(:,3)*100, 'color','blue', 'linewidth',4)
hold on
zzz = squeeze(dsge_irf(3,:,4));
plot(0:hMax,zzz*100, 'linewidth',4, 'color','black','linestyle','--')
zzz = squeeze(dsge_irf_new(3,:,4));
plot(0:hMax,zzz*100, 'linewidth',4, 'color','red','linestyle',':')
ylim([-.9, .1])
xlim([0,6])
xlabel('Horizon h','fontsize',22)
set(gca,'fontsize',20,'fontweight','demi')
box off
hold off
print('-dpng',[picPath, 'figure2_3'])

%dashed (black) line is DSGE IRF (zeta = .65), solid (blue) is VAR implied,
%dotted (red) lines is DSGE (zeta = .5)


%% Figure 2
%%%%  Match data VAR irf to DSGE IRF

numReps = 100;                           % number of simulations
sampleSizes = [80, 200];                 % sample sizes to consider
Tau = linspace(.2,.95, 50);                % vary zetaP across the grid
bbeta = theta(1);                        % parameter to calculate monetary policy rul
zetaP = theta(5);                        % to plot the true zetaP
hMax  = 10;
zeta_opt = zeros(numReps,length(sampleSizes));     %store all optimal zetas      
zeta_opt2 = zeta_opt;
err = zeros(2,10);


parfor n = 1:length(sampleSizes)
    
    N = sampleSizes(n);
    rng(n*100)

    for jj = 1:numReps
        
        [Yhat, s_true] = estimDSGE(theta, N, 100);  %simulate data
        
        %calculate VAR matrices
        Y_true = [Yhat(:,4) - Yhat(:,3)/bbeta, Yhat(:,1), Yhat(:,3)];
        p = 1;
        
        Yp = Y_true(1:length(Y_true)-p,:);
        Xp = Y_true(2:length(Y_true)-(p-1),:);
        Xp = [ones(length(Xp),1),Xp];
        
        PhiHat_true = inv(Xp'*Xp)*Xp'*Yp;
        sigmaHat = 1/length(Yp)*(Yp - Xp*PhiHat_true)'*(Yp - Xp*PhiHat_true);
        Sigma_chol = chol(sigmaHat)';
        
        F = PhiHat_true(2:length(PhiHat_true),:)';
        VAR_irf =  zeros(hMax+1,3);
        for h = 0:hMax
            VAR_irf(h+1,:) = F^h * Sigma_chol(:,1);
        end
        
        Q = zeros(length(Tau),1);
        
        theta_new = theta;
        err = zeros(2, size(VAR_irf,1));
        for z = 1:length(Tau)
            zeta = Tau(z);
            theta_new(5) = zeta;
            
            %calculate DSGE implied IRFs
            irf = zeros(4,hMax+1,5);
            for j = 1:4
                for h = 0:hMax
                    irf(:,h+1,j) = ImpulseResponse(theta_new,h,j);
                end
            end
            
            j = 4;
            zzz = squeeze(irf([1,3],:,j));

            err(1,:) = cumsum(VAR_irf(:,2)) - cumsum(zzz(1,:))';
            err(2,:) = VAR_irf(:,3) - zzz(2,:)';
            
            Q(z) = mean( mean(err.^2) );
            
        end
        
        [~,minLoc] = min(Q);
        zeta_opt(jj,n) = Tau(minLoc);
        
        %match on DSGE representation of VAR 
        theta_new = theta;
        for z = 1:length(Tau)
            zeta = Tau(z);
            theta_new(5) = zeta;
            
            %calculate VAR approximation to DSGE IRFs
            
            [Phi_star, Sigma_star] = DSGE_VAR_3(theta_new, 1);
            
            Sigma_chol = chol(Sigma_star)';
            F = Phi_star(1:length(Phi_star)-1,:)';   %companion Form
            
            m = size(Sigma_chol,1);

            M_n =[eye(m), zeros(m,(p-1)*m)];
            C_hat = @(h) M_n*F^h*M_n';
            
            
            irf = zeros(hMax+1,m);
            for h = 0:hMax
                irf(h+1,:) = C_hat(h) * Sigma_chol(:,1);
            end
            
            err(1,:) = cumsum(VAR_irf(:,2)) - cumsum(irf(:,2));
            err(2,:) = VAR_irf(:,3) - irf(:,3);
            
            Q(z) = mean( mean(err.^2) );
            
        end
        
        [~,minLoc] = min(Q);     
        zeta_opt_2(jj,n) = Tau(minLoc);
    end
end

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

print('-dpng',[picPath, 'DSGE_irf_matching'])

hold off


[d, grid] = ksdensity(zeta_opt_2(:,1));
plot(grid,d, 'linewidth',4, 'Color','blue', 'linestyle',':')
hold on
[d, grid] = ksdensity(zeta_opt_2(:,2));
plot(grid,d, 'linewidth',4, 'Color','black', 'linestyle','--')

plot([zetaP zetaP], [0 max(d)+10],...
        'Color','red','LineWidth',2.5 ), hold off
ylim([0,8])
xlim([0,1])
set(gca, 'Xtick',0:.2:1)
set(gca,'fontsize',20,'fontweight','demi')

print('-dpng',[picPath, 'DSGE_VAR_matching'])

hold off

%red (solid) is true value. Black (dashed) is sample size 200, while blue
%(dotted) is sample size 80