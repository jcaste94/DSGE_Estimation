%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                        DSGE ESTIMATION                    
%
% Author: Juan Castellanos Silván building on Jacob Warren's code
% Date  : 24/04/2020
%==========================================================================

% Problem 2: Filtered states (Fig 24 in FVRRS Handbook Chapter)

%==========================================================================
%                          HOUSEKEEPING
%==========================================================================

tic
clear all
close all
clc

rng(123);
% add helper functions:
% addpath(strcat(pwd,'/Replication Code (MATLAB)/helper_functions'));

%==========================================================================
%                   GENERATE DATA FROM SOLVED DSGE
%==========================================================================

% pass in parameters from paper
theta = parameters;

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);

% generate time series with 50 observations (100 used as burn-in to get to
% steady state
[Yhat, s_true] = estimDSGE(theta, 50, 100);


%==========================================================================
%                           RUN KALMAN FILTER
%==========================================================================

% -----------------
% 0. Initialization
% ------------------
Sigma_u = zeros(4);
alpha = .1;

% ------------------------------------
% 1. Run the filter based on only Y(1)
% -------------------------------------
[s_mean1, P_mean1, s_hi1, s_lo1, L1] = kal_wrapper(theta,Yhat(:,1),alpha, Sigma_u);


% ------------------------------
% 2. Run the filter on Y(:,1:3)
% ------------------------------
[s_mean13, P_mean13, s_hi13, s_lo13, L13] = kal_wrapper(theta,Yhat(:,1:3),alpha, Sigma_u);


%==========================================================================
%                               FIGURES
%==========================================================================

% helper values to print the background
x = 1:50;              
X = [x,fliplr(x)];

state_list = {'preference', 'markup', 'technology', 'monetary','lag_output'};

% Only output Y(1) 
 for i = [1 3]   
    
    figure('Name',state_list{i})
    hold on
    Y_fill = [s_lo1(:,i)',fliplr(s_hi1(:,i)')];
    h = fill(X,Y_fill,[0 0 0] + .8);
    set(h, 'EdgeColor','None')
    plot(s_true(:,i), 'LineWidth',4,'Color','k', 'LineStyle',':')
    plot(s_mean1(:,i),'Color','b','LineStyle','--','LineWidth',4)
    set(gca,'fontsize',20,'fontweight','demi')
    grid on 
    xlabel('Time','fontsize',22)
    xlim([0,50])
    if i==1
        ylim([-.06 .1])
    else
        ylim([-.03 .04])
    end
    
    box off
    
    x = 29.7;                  % A4 paper size
    y = 21.0;                  % A4 paper size
    xMargin = 1;               % left/right margins from page borders
    yMargin = 1;               % bottom/top margins from page borders
    xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
    ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

    set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

    set(gcf, 'PaperUnits','centimeters')
    set(gcf, 'PaperSize',[x y])
    set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
    set(gcf, 'PaperOrientation','portrait')
    
    cd('/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS3/LaTeX/')
    saveas(gcf, strcat('pFS1_', state_list{i},'.pdf'));

 
 end
 
 
 % Add inflation and labor share Y(1:3)
 for i = [1 3]   
    
    figure('Name',state_list{i})
    hold on
    Y_fill = [s_lo13(:,i)',fliplr(s_hi13(:,i)')];
    h = fill(X,Y_fill,[0 0 0] + .8);
    set(h, 'EdgeColor','None')
    plot(s_true(:,i), 'LineWidth',4,'Color','k', 'LineStyle',':')
    plot(s_mean13(:,i),'Color','b','LineStyle','--','LineWidth',4)
    set(gca,'fontsize',20,'fontweight','demi')
    grid on 
    xlabel('Time','fontsize',22)
    xlim([0,50])
    if i==1
        ylim([-.06 .1])
    else
        ylim([-.03 .04])
    end
    
        box off
    
    x = 29.7;                  % A4 paper size
    y = 21.0;                  % A4 paper size
    xMargin = 1;               % left/right margins from page borders
    yMargin = 1;               % bottom/top margins from page borders
    xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
    ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

    set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

    set(gcf, 'PaperUnits','centimeters')
    set(gcf, 'PaperSize',[x y])
    set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
    set(gcf, 'PaperOrientation','portrait')
    
    cd('/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS3/LaTeX/')
    saveas(gcf, strcat('pFS13_', state_list{i},'.pdf'));

 end
 
 cd('/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS3/DSGE_estimation/')

toc