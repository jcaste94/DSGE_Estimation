%--------------------------------------------
%--------------------------------------------
% File to perform DSGE implied kalman and particle filtering
%
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
picPath = [pwd, '/figuresv2/likelihood/'];

% pass in parameters from paper
theta = parameters;

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);



%% Kalman Filtering

%generate time series with 50 observations (100 used as burn-in to get to
%steady state
[Yhat, s_true] = estimDSGE(theta, 50, 100);

%run the filter/smoother on Y(:,1:3)
Sigma_u = zeros(4);
alpha = .1;
[s_mean, P_mean, s_hi, s_lo, L] = kal_wrapper(theta,Yhat(:,1:3),alpha, Sigma_u);

%plot the true state and the filtered state
x = 1:50;                %helper values to print the background
X = [x,fliplr(x)];
for i = [1 3]
    figure(i)
    hold on
    Y_fill = [s_lo(:,i)',fliplr(s_hi(:,i)')];
    h = fill(X,Y_fill,[0 0 0] + .8);
    set(h, 'EdgeColor','None')
    plot(s_true(:,i), 'LineWidth',4,'Color','black', 'LineStyle',':')
    plot(s_mean(:,i),'Color','b','LineStyle','--','LineWidth',4)
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Time','fontsize',22)
    
    xlim([0 50])
    if i==1
        ylim([-.06 .1])
    else
        ylim([-.03 .04])
    end
    box off
    print('-dpng',[picPath, 'filtered' num2str(i) 'from_3'])
    
    hold off
end

close all


%Now filter based on only Y(1)
[s_mean, P_mean, s_hi, s_lo, L] = kal_wrapper(theta,Yhat(:,1),alpha, Sigma_u);

%plot the true state and the filtered state
for i = [1 3]
    figure(i)
    hold on
    Y_fill = [s_lo(:,i)',fliplr(s_hi(:,i)')];
    h = fill(X,Y_fill,[0 0 0] + .8);
    set(h, 'EdgeColor','None')
    plot(s_true(:,i), 'LineWidth',4,'Color','black', 'LineStyle',':')
    plot(s_mean(:,i),'Color','b','LineStyle','--','LineWidth',4)
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Time','fontsize',22)
    
    xlim([0 50])
    if i==1
        ylim([-.06 .1])
    else
        ylim([-.03 .04])
    end
    box off
    print('-dpng',[picPath, 'filtered' num2str(i) 'from_1'])
    
    hold off
end

close all

%black dotted line is true state, blue dashed line is kalman filtered
%state. Grey band represents 90% credible bands

%% Particle Filter

Sigma_u = diag(var(Yhat(:,1:3)) * .1);


%first, re-estimate Kalman Filter with all observables
alpha = .1;
[s_mean, P_mean, s_hi, s_lo,  L] = kal_wrapper(theta,Yhat(:,1:3),alpha, Sigma_u);

%Need estimation error now in order to evaluate the particles
resample = 1;                %whether to do resampling in the particle filter
N_options = [100, 200, 500];  %Number of particles
Nsim = 100;                    %Number of simulations to compare with
s_i = zeros(Nsim,50,5);
L_i = zeros(Nsim,1);


for n = 1:length(N_options)
    close all
    
    N = N_options(n)
    
    %Plot density of particle filter likelihood approximation minus kalman
    %filter true likelihood
    for i = 1:Nsim
        [s_ret,L_ret,N_eff_ret,hi,lo] = particle_wrapper(theta,Yhat(:,1:3),N,resample,Sigma_u,alpha);
        L_i(i) = sum(L_ret(2:end));
    end
    [f,xi] = ksdensity(L_i -sum(log(L(2:end))) );
    
    f_part(n,:) = f;    %save for plotting all densities together
    xi_part(n,:) = xi;  %save for plotting all densities together
    
    
    %only do plots for N = 100
    if N==100
        
        for i = 1:4
            plot(s_mean(:,i), 'color','red','Linewidth',4)
            hold on
            
            plot(s_ret(:,i), 'color','blue','Linewidth',4,'linestyle','--')
            
            set(gca,'fontsize',20,'fontweight','demi')
            xlabel('Time','fontsize',22)
            
            xlim([0 50])
            
            hold off
            box off
            
            print('-dpng',[picPath, 'Particle_filtered' num2str(N) '_' num2str(i) ])
        end
        
        
        figure(2)
        hold on
        
        plot(log(L), 'color','red','Linewidth',4)
        plot(L_ret, 'color','blue','Linewidth',4,'linestyle','--')
        hold off
        set(gca,'fontsize',20,'fontweight','demi')
        xlabel('Time','fontsize',22)
        xlim([0 50])
        
        print('-dpng',[picPath 'Particle_likelihood' num2str(N) ])
        
    end
end

close all

p = plot(xi_part',f_part','Linewidth',4);
p(2).LineStyle = ':';
p(3).LineStyle = '--';
box off
set(gca,'fontsize',20,'fontweight','demi')
xlabel('Error','fontsize',22)
%set(gca,'XTick', -4:1:2 )
print('-dpng',[picPath 'Particle_L_together'])

%blue solid is 100, red dotted is 200, yellow dashed is 500
