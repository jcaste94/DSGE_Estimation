%--------------------------------------------
%--------------------------------------------
% File to perform DSGE implied autocovariances, spectrum, forecast error
% variance decompositions, impulse response functions,
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
picPath = [pwd, '/figuresv2/modelimplications/'];

% pass in parameters from paper
theta = parameters;

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);


%% Autocovariances
h = 0;
[Gamma_ss_h, Gamma_yy_h, Gamma_ys_h] = DSGEautocovar(theta, h);

yy_corr = corrcov(Gamma_yy_h);

H = 12;

autocorr_plot = zeros(4,H+1);
for h = 0:H
    [~, yy_cov, ~] = DSGEautocovar(theta, h);
    autocorr_plot(1,h+1) = yy_cov(1,1)/Gamma_yy_h(1,1);
end

H2 = linspace(-6,6,13);
for h = 1:length(H2)
    [~, yy_cov, ~] = DSGEautocovar(theta, H2(h));
    
    autocorr_plot(2,h) = yy_cov(1,2)/sqrt(Gamma_yy_h(1,1)*Gamma_yy_h(2,2));
    autocorr_plot(3,h) = yy_cov(1,3)/sqrt(Gamma_yy_h(1,1)*Gamma_yy_h(3,3));
    autocorr_plot(4,h) = yy_cov(1,4)/sqrt(Gamma_yy_h(1,1)*Gamma_yy_h(4,4));
    
end

%plot output autocorrelation
plot(0:H,autocorr_plot(1,1:13), 'linewidth',4, 'Color','blue')
xlim([0,12])
ax = gca;
set(ax,'XTick',[0:2:12])
set(gca,'fontsize',20,'fontweight','demi')
xlabel('Temporal Shift h','fontsize',22)
ylim([-.3 1.2])
box off

print('-dpng',[picPath, 'autocorr_' num2str(1)])

%plot the three ACF functions together:
close all
figure(1)
hold on
plot(-6:6,autocorr_plot(2,:),'linewidth',4, 'linestyle','-','Color','blue')
plot(-6:6,autocorr_plot(3,:),'linewidth',4, 'linestyle',':','Color','black')
plot(-6:6,autocorr_plot(4,:),'linewidth',4, 'linestyle','--','Color',[0 .8 0])
xlim([-6,6])
ax = gca;
set(ax,'XTick',[-6:2:6])
hold off
set(gca,'fontsize',20,'fontweight','demi')
xlabel('Temporal Shift h','fontsize',22)
print('-dpng',[picPath, 'autocorr_with_output'])

%first is labor, second inflation, third interest rate


%% Forecast error Variance decomposition

H = linspace(1, 20, 20);

Y_err = zeros(4,4,length(H));

for h = 1:length(H)
    [FEVD, ~] = DSGEforecasterror( theta,H(h) );
    Y_err(:,:,h) = FEVD;
end

for i = 1:4
    zzz = squeeze(Y_err(i,:,:));    
    
    b= bar(zzz','stack');
    axis([1 max(H) 0 1])
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Horizon h','fontsize',22)
    
    C= [ [0,0,0]; [0,0,0]+.3; [0,0,0] + .6; [0,0,0]+.9 ];
    
    for n = 1:length(b)
        b(n).FaceColor = C(n,:);
    end
    
    ax = gca;
    set(ax,'XTick', 0:5:max(H) )
    
    print('-dpng',[picPath, 'FEVD_' num2str(i)])
end

%darkest is phi_t follows by lambda_t, z_t, eps_R_t.

%% Spectrum

omega = linspace(.196, .785, 50);

[~, ~, ~, f_yy] = DSGEspectrum(theta, omega);


for i = 1:4
    zzz = squeeze(f_yy(i,1:4,:));
        
    b = bar(omega,zzz','stack', 'EdgeColor','black');
    xlim([min(omega) max(omega)])
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Frequency \omega','fontsize',22)
    
    
    C= [ [0,0,0]; [0,0,0]+.3; [0,0,0] + .6; [0,0,0]+.9 ];
    
    for n = 1:length(b)
        b(n).FaceColor = C(n,:);
        b(n).EdgeColor = 'black';
    end
    
    box off
    
    ax = gca;
    set(ax,'XTick', .2:.1:max(omega) )
    
    
    print('-dpng',[picPath, 'spectrum_stacked' num2str(i)])
end

%black is shock1 through shock 4. Note shock 5 is not plotted because it is
%linear combination of other shocks


%% Impulse Response Functions

hMax = 20;
irf = zeros(4,hMax+1,5);
for j = 1:4
    for h = 0:hMax
        irf(:,h+1,j) = ImpulseResponse(theta,h,j);
    end
end

%plot IRF to log output - order is 1-preferences, 2-mark-up, 3-technology,
%4-monetary policy
for j = 1:4
    zzz = squeeze(irf(1,:,j));
    plot(0:hMax,cumsum(zzz*100), 'linewidth',4, 'color','blue')
    
    if j==4
        ylim([-.85 .2])
    elseif j==1 || j==2
        ylim([-.8 0.1])
    elseif j==3
        ylim([-.1 1.2])
    end
    if j== 1 || j==2
        hold on
        plot(0:hMax, zeros(hMax+1,1), 'linewidth', 2.5, 'color',[0 0 0 ]+.5, ...
            'linestyle',':')
        hold off
    end
    box off
    
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Horizon h','fontsize',22)
    print('-dpng',[picPath, 'IRF_' num2str(1) '_' num2str(j)])
end
