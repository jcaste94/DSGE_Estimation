function [L,s_t,postWeights,N_eff] = particle_filter(Y,Psi0,...
    Psi1, Phi_eps, Phi1,Sigma_u,particles,weights,resample)

%Perform filtering for date t

T = length(Y);
[N, p] = size(particles);

%draw one-step ahead forecast particles
eps_t = mvnrnd(zeros(p,1),Phi_eps*Phi_eps',N);

s_tilde = particles*Phi1' + eps_t;

y_hat = ones(N,1)*Psi0 + s_tilde*Psi1';

v = ones(N,1)*Y-y_hat;     %forecast error

%define un-normalized weights
weights_tilde =weights.*mvnpdf(v,0,Sigma_u); 

%approximate log-likelihood
L = log(mean(weights_tilde)); 

%define normalized weights
pi_t = weights_tilde./(mean(weights_tilde));

%effective number of particles
N_eff = N^2/sum(pi_t.^2);

%resample randomly if the effective sample size is too small
%(conditionally)
s_t = s_tilde;
postWeights = pi_t;

if resample==1 && N_eff<.5*N
    randInts = randsample(1:N,N,true,pi_t);
    s_t = s_tilde(randInts,:);
    postWeights=ones(N,1);
end


