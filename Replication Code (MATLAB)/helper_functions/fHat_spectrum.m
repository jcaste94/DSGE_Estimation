function [f,f_bar] = fHat_spectrum(Y, omega, a, c)

[T,n] = size(Y);


autocovs = zeros(2*T-1,n);
for i = 1:n
    autocovs(:,i) = xcov(Y(:,i),'biased');
end


%H = ( -(T-1):T-1 )';
H = (0:T-1)';
j = -(T-1)/4:(T-1)/2;   %1:(T-1)/2;

%K = @(omega_j,omega,lambda) abs((omega_j-omega)/lambda)<1/2;

K = @(omega_j, omega, a,c) a*exp(-(omega-omega_j).^2 ./ (2*c^2));

omega_j = j*2*pi/T;
f = zeros(length(omega_j),n);
autocovs = autocovs(T:end,:);
gamma0 = autocovs(1);
gamma_h = autocovs(2:end,:);
for o = 1:length(omega_j)
    for i = 1:n
        cos_term = cos(omega_j(o).*H);
        f(o,i) = 1/(2*pi)*(sum(2*autocovs(:,i).*(cos_term) ) - autocovs(1,i)) ;
    end
   %f(o,:) = 1/(2*pi)*(gamma0 + 2*sum(gamma_h.*cos(omega(o).*H)));
end 

f_bar = zeros(length(omega),n);
for o = 1:length(omega)
   %f_bar(o,:) = pi/(lambda*(T-1)/2) * K(omega_j,omega(o),lambda)*f ;
   f_bar(o,:) = pi/(c*(T-1)/2) * K(omega_j,omega(o),a, c)*f ;
end 
   

end