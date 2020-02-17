
%--------------------------------------------
%--------------------------------------------
% File to estimate empirical VAR, autocovariances, spectrum, and impulse
% response functions
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
picPath = [pwd, '/figuresv2/empiricalanalogues/'];

% pass in parameters from paper
theta = parameters;

%% Estimate VAR by data

load Data/Hours
load Data/FedFunds
load Data/GDP
load Data/GDPdef
load Data/Labor

Y = [GDPdat LabShare GDPDEF FedFundsQuarterly];

%% Use BIC to find optimal lag length

maxLag = 6;
for i = 1:maxLag
    Yp = Y(1:length(Y)-i,:);
    Xp = Y(2:length(Y)-(i-1),:);
    for j = 2:i
        Yp = [Yp, Y(j:length(Y)-(i-j+1),:)];
        Xp = [Xp, Y(j+1:length(Y)-(i-j),:)];
    end
    
    Xp = [ones(length(Xp),1),Xp];
    
    PhiHat = inv(Xp'*Xp)*Xp'*Yp;
    sigmaHat = 1/length(Yp)*(Yp - Xp*PhiHat)'*(Yp - Xp*PhiHat);
    sigmaHat = sigmaHat(1:4,1:4);
    
    
    BIC(i) = length(Yp)*log(det(sigmaHat)) + (1+size(Y,2)^2*i)/2*log(length(Yp));
end

[c,lagLength] = min(BIC);
%optimal lag length is 1 by BIC

i = lagLength;
Yp = Y(1:length(Y)-i,:);
Xp = Y(2:length(Y)-(i-1),:);
for j = 2:i
    Yp = [Yp, Y(j:length(Y)-(i-j+1),:)];
    Xp = [Xp, Y(j+1:length(Y)-(i-j),:)];
end

Xp = [ones(length(Xp),1),Xp];

PhiHat = inv(Xp'*Xp)*Xp'*Yp;
F = PhiHat(2:length(PhiHat),:)';

sigmaHat = 1/length(Yp)*(Yp - Xp*PhiHat)'*(Yp - Xp*PhiHat);

Psi0Hat = PhiHat(1,1:4);
Psi1Hat = PhiHat(2:length(PhiHat),1:4)';


SigmaY = 1/length(Yp)*Yp'*Yp;

vecSigma = inv(eye(size(SigmaY,2)^2) - kron(F,F))*sigmaHat(:);
Sigma0 = reshape(vecSigma,size(SigmaY,2),size(SigmaY,2));

M_n = [eye(4), zeros(4,4*lagLength-4)];
Gamma0 = M_n*Sigma0*M_n';

%% Autocovariance Function

Gamma_j = zeros(4,4,25);
count = 0;
for i = -12:12
    count = count+1;
    if(i>0)
        Gamma_j(:,:,count) = M_n*(F^i)*Sigma0*M_n';
    elseif i<0
        Gamma_j(:,:,count) = (M_n*((F^abs(i))*Sigma0)*M_n')';
    else
        Gamma_j(:,:,count) = (M_n*(Sigma0)*M_n');
    end
end

for i = 13:25
    autocorr_plot(1,i-12) = Gamma_j(1,1,i)/Gamma0(1,1);
end

for h = 7:19
    autocorr_plot(2,h-6) = Gamma_j(1,2,h)/sqrt(Gamma0(1,1)*Gamma0(2,2));
    autocorr_plot(3,h-6) = Gamma_j(1,3,h)/sqrt(Gamma0(1,1)*Gamma0(3,3));
    autocorr_plot(4,h-6) = Gamma_j(1,4,h)/sqrt(Gamma0(1,1)*Gamma0(4,4));
    
end

close all
hold on
plot(-6:6,autocorr_plot(2,:),'linewidth',4, 'linestyle',':','Color','blue')
plot(-6:6,autocorr_plot(3,:),'linewidth',4, 'linestyle','--','Color','red')
plot(-6:6,autocorr_plot(4,:),'linewidth',4,'Color','black')
hold off
xlim([-6,6])
set(gca, 'Xtick',-6:2:6)
ylim([-.4,.4])
set(gca, 'Ytick', -.4:.2:.4)
set(gca,'fontsize',20,'fontweight','demi')
xlabel('Temporal Shift h','fontsize',22)
print('-dpng',[picPath, 'Empirical_autocorr_with_output'])
%black (solid) is interest rate, red (dashed) is inflation, blue (dotted)
%is labor


%sample autocorrelation
full_sample = zeros(13,3);
for i = 2:4
    
    xcf = crosscorr(GDPdat,Y(:,i), 6);
    
    full_sample(:,i-1) = xcf;
end

close all
hold on
plot(-6:6,full_sample(:,1),'linewidth',4, 'linestyle',':','Color','blue')
plot(-6:6,full_sample(:,2),'linewidth',4, 'linestyle','--','Color','red')
plot(-6:6,full_sample(:,3),'linewidth',4,'Color','black')
hold off
xlim([-6,6])
set(gca, 'Xtick',-6:2:6)
set(gca, 'Ytick', -.4:.2:.4)
set(gca,'fontsize',20,'fontweight','demi')
xlabel('Temporal Shift h','fontsize',22)
print('-dpng',[picPath, 'sample_autocorr_with_output'])

%% Spectrum


%plot figure by figure, VAR implied and then smoothed periodogram. Figured
%are in order of 1-output, 2-labor, 3- inflation 4-interest rates.
%I use stdev = .15 for gaussian kernel. Dotted line is VAR implied, solid
%is smoothed periodogram

M = @(z) [eye(4)*z];

i = sqrt(-1);
omega = linspace(0, 1, 500);
omega_BS = linspace(.196, .785, 50);

f_yy = zeros(4,length(omega));
for o = 1:length(omega)
    Sl = eye(4) - Psi1Hat*M(exp(i*omega(o)))';
    f_yy(:,o) = real(diag(1/(2*pi)*inv(Sl)*sigmaHat*inv(Sl')));
end

[f,f_bar] = fHat_spectrum(Y,omega, 1, .15);
f_bar(f_bar<=0) = .0001;

for i = 1:4
    hold off
    plot(omega,f_bar(:,i), 'linewidth',4, 'Color','blue');
    hold on
    
    max_f = max(f_bar(:,i)) + sqrt(var(f_bar(:,i)));
    patch([min(omega_BS) max(omega_BS) max(omega_BS) min(omega_BS)], [0,0,max_f,max_f], .9*[1,1,1])
    plot(omega,f_bar(:,i), 'linewidth',4, 'Color','blue')
    plot(omega,f_yy(i,:), 'linewidth',4, 'Color','red', 'LineStyle',':');
    ylim([0, max_f])
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Frequency \omega','fontsize',22)
    box off
    
    print('-dpng',[picPath, 'Spectrum' num2str(i)])
    
end


%% Impulse Responses
% sign restrictions: inflation must be negative and interest rate must be
% positive for the first four quarters

[m,n] = size(PhiHat);
hMax = 30;
irf = zeros(hMax+1,4);

k=1;
for i = 1:1000
    e = randn(4,1);
    
    Omega = (e./norm(e));
    
    Sigma_tr = chol(Gamma0)';
    C_hat = @(h) M_n*F^h*M_n';
    
    for j = 0:hMax
        irf(j+1,:) = (C_hat(j)*Sigma_tr*Omega)';
    end
    
    if irf(1:4,3)<=0 & irf(1:4,4)>=0
        irf_keep(k,:,:) = irf;
        Omega_keep(k,:) = Omega;
        k = k+1;
    end
    
end

irfMin = zeros(hMax+1,n);
irfMax = zeros(hMax+1,n);
for i = 1:hMax+1
    irfMin(i,:) = min(squeeze(irf_keep(:,i,:)));
    irfMax(i,:) = max(squeeze(irf_keep(:,i,:)));
end

%find specific irf such that output falls initialls
q_select = find(squeeze(irf_keep(:,1,1)) < 0, 1);


Omega_low = min(Omega_keep);
Omega_high = max(Omega_keep);


namesY = {'GDP','LabShare','GDPDEF','FedFunds'};
x = 0:hMax;
X = [x,fliplr(x)];


for i = 1:4
    hold off
    if(strcmp(namesY{i},'GDP'))
        
        Y_fill = [cumsum(100*irfMin(:,i))',fliplr(cumsum(100*irfMax(:,i))')];
        h = fill(X,Y_fill,[0 0 0] + .8);
        set(h, 'EdgeColor','None')
        hold on
        plot(0:hMax, cumsum(100*irf_keep(q_select,:,i)), 'linewidth',4,'color','blue')
        
        
        %         plot(0:hMax, cumsum(100*irf05(:,i)), 'linewidth',4, 'color','blue')
        %         hold on
        %         plot(0:hMax, cumsum(100*irf95(:,i)), 'linewidth',4, 'color','blue')
        %
    else
        
        Y_fill = [100*irfMin(:,i)',fliplr(100*irfMax(:,i)')];
        h = fill(X,Y_fill,[0 0 0] + .8);
        set(h, 'EdgeColor','None')
        hold on
        plot(0:hMax, 100*irf_keep(q_select,:,i), 'linewidth',4,'color','blue')
        %
        %
        %         plot(0:hMax, 100*irf05(:,i), 'linewidth',4, 'color','blue')
        %         hold on
        %         plot(0:hMax, 100*irf95(:,i), 'linewidth',4, 'color','blue')
    end
    
    set(gca,'fontsize',20,'fontweight','demi')
    xlabel('Horizon h','fontsize',22)
    box off
    
    print('-dpng',[picPath, 'irf_signrestrictions', num2str(i)])
    
end



