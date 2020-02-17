function parasim = priordraws(nsimul)
% Function to generate prior draws
% - Now prior draws are within bounds
% Minchul Shin (Penn)
% Last Date: 5/8/2014



% The priors considered are:
%  1.) 100*(1/Beta-1) is GAMMA with mean 0.5 and std 0.5
%  2.) 100*log(Gamma) is NORMAL with mean .75 and std 0.5
%  3.) Lambda is GAMMA with mean .2 and std 0.2
%  4.) 100*log(PiStar) is GAMMA with mean 1 and std 0.5
%  5.) zetaP is BETA with mean 0.7 and std 0.15
%  6.) 1/(1+nu) is GAMMA with mean 1.50 and std 0.75
%  7.) rhoPhi is UNIFORM with mean 0 and std 1
%  8.) rhoLambda is UNIFORM with mean 0 and std 1
%  9.) rhoZ is UNIFORM with mean 0 and std 1
%  10.) 100*sigmaPhi is InvGamma with mean 2 and std 4.0
%  11.) 100*sigmaLambda is InvGamma with mean 0.5 and std 4.0
%  12.) 100*sigmaZ is InvGamma with mean 2 and std 4.0
%  13.) 100*sigmaR is InvGamma with mean 0.5 and std 4


% define some variables
npara = 13;

parasim = zeros(nsimul,npara);


% % GAMMA pdf
para1 = [0.5,  0.2, 1, 1.50];
para2 = [0.5,  0.2,  0.5, 0.75];

b = para2.^2./para1;
a = para1./b;

sims = zeros(nsimul,4);
for i = 1:length(para1)
    b = para2(i)^2/para1(i);
    a = para1(i)/b;
    
    sims(:,i) = gamrnd(a,b,nsimul,1);
end

parasim(:,[1,3,4, 6]) = sims;


% % NORMAL pdf
para1 = 0.75;
para2 = 0.5;

parasim(:,2) = para1 + para2*randn(nsimul,1);


% % BETA pdf
para1 = 0.7;
para2 = 0.15;

a = (1-para1).*para1.^2./para2.^2 - para1;
b = a.*(1./para1 - 1);

parasim(:,5) = betarnd(a,b,nsimul,1);


% % UNIFORM pdf

parasim(:,7:9) = unifrnd(0,1,[nsimul,3]);


%Inverse gamma
para1 = [2, 0.5, 2, 0.5]';
para2 = [4, 4, 4, 4]';

a = para1;
b = para2;

for i = 1:length(para1)
    b = para2(i);
    a = para1(i);
    parasim(:,9+i) = sqrt( b*a^2./sum( (randn(b,nsimul)).^2 )' );
end

            
% 
% 
% 
% 
% 
% % some housekeeping
% pshape   = prior(:,1);
% pmean    = prior(:,2);
% pstdd    = prior(:,3);
% pmask    = prior(:,4);
% pfix     = prior(:,5);
% pmaskinv = 1- pmask;
% pshape   = pshape.*pmaskinv;
% 
% % parameter specification: transformation
% trspec(:,1) = trspec(:,1).*pmaskinv;
% bounds      = trspec(:,2:3);
% 
% % define some variables
% npara = size(prior,1);
% 
% % loop to generate parameter draws
% parasim = zeros(nsimul,npara);
% 
% for i=1:npara
%     
%     if pmask(i)==1
%         parasim(:,i) = ones(nsimul,1)*pfix(i);
%         
%     else
%         
%         % beta prior
%         if pshape(i)==1
%             a = (1-pmean(i))*pmean(i)^2/pstdd(i)^2 - pmean(i);
%             b = a*(1/pmean(i) - 1);
%             parasim(:,i) = betarnd(a,b,nsimul,1);
%             
%             % check bound and re-draw
%             for j=1:1:sum(outbound)
%                 notvalid = 1;
%                 while notvalid
%                     temp_para = betarnd(a,b,1,1);
%                     
%                     if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
%                         parasim(outboundset(j),i) = temp_para;
%                         notvalid = 0;
%                     end
%                 end
%             end
%             
%             % gamma prior
%         elseif pshape(i) == 2
%             b = pstdd(i)^2/pmean(i);
%             a = pmean(i)/b;
%             parasim(:,i) = gamrnd(a,b,nsimul,1);
%             
%             % check bound and re-draw
%             for j=1:1:sum(outbound)
%                 notvalid = 1;
%                 while notvalid
%                     temp_para = gamrnd(a,b,1,1);
%                     
%                     if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
%                         parasim(outboundset(j),i) = temp_para;
%                         notvalid = 0;
%                     end
%                 end
%             end
%             
%             % gauusian prior
%         elseif pshape(i)==3
%             a = pmean(i);
%             b = pstdd(i);
%             parasim(:,i) = a + b*randn(nsimul,1);
%             
%             % check bound and re-draw
%             for j=1:1:sum(outbound)
%                 notvalid = 1;
%                 while notvalid
%                     temp_para = a + b*randn(1,1);
%                     
%                     if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
%                         parasim(outboundset(j),i) = temp_para;
%                         notvalid = 0;
%                     end
%                 end
%             end
%             
%             % inverse gamma prior
%         elseif pshape(i)==4
%             a = pmean(i);
%             b = pstdd(i);
%             parasim(:,i) = sqrt( b*a^2./sum( (randn(b,nsimul)).^2 )' );
%             
%             % check bound and re-draw
%             outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
%             outboundset = find(outbound);
%             for j=1:1:sum(outbound)
%                 notvalid = 1;
%                 while notvalid
%                     temp_para = sqrt( b*a^2./sum( (randn(b,1)).^2 )' );
%                     
%                     if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
%                         parasim(outboundset(j),i) = temp_para;
%                         notvalid = 0;
%                     end
%                 end
%             end
%             
%             % uniform prior
%         elseif pshape(i)==5
%             a = pmean(i);
%             b = pstdd(i);
%             parasim(:,i) = a + (b-a)*rand(nsimul,1);
%             
%             % check bound and re-draw
%             outbound = ~((bounds(i,1) < parasim(:,i)).*(parasim(:,i) < bounds(i,2)));
%             outboundset = find(outbound);
%             for j=1:1:sum(outbound)
%                 notvalid = 1;
%                 while notvalid
%                     temp_para = a + (b-a)*rand(1,1);
%                     
%                     if ((bounds(i,1) < temp_para)&&(temp_para < bounds(i,2)))
%                         parasim(outboundset(j),i) = temp_para;
%                         notvalid = 0;
%                     end
%                 end
%             end
%             
%             % no prior, fixed
%         elseif pshape(i)==0
%             a = pmean(i);
%             parasim(:,i) = a*ones(nsimul,1);
%             
%         end
%     end
% end	% i loop
% 












