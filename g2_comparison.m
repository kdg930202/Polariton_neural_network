function [g2_s_an, g2_s_nu] = g2_comparison(phi, theta)

p1 = cos(theta)^2;
p2 = sin(theta)^2 * cos(phi)^2;
p3 = sin(theta)^2 * sin(phi)^2;

d = 8; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
p_number = 1;
I = eye(d);
alpha = 1;

%%
% Coherent state
rho_coh = 0; %Initialization

for i = 0:d-1
    for j = 0:d-1
        rho_coh = rho_coh + exp(-alpha^2)*(alpha^i)*(conj(alpha)^j)/sqrt(prod(1:i)*prod(1:j))*I(:,i+1)*I(:,j+1)';
    end
end

% trace(rho_coh)
g2_coh = trace(a'*a'*a*a*rho_coh)/trace(a'*a*rho_coh)^2; %should be 2

%% 
% Thermal state
nth = 0.5;
rho_th = 0;

for i=0:d-1
    rho_th = rho_th+(nth^i/(1+nth)^(i+1))*I(:,i+1)*I(:,i+1)'; 
end
% trace(rho_th)
g2_th = trace(a'*a'*a*a*rho_th)/trace(a'*a*rho_th)^2; %should be 1

%%
% rho_fock = I(:,d)*I(:,d)';
psi_fock = (a')^p_number/sqrt(prod(1:p_number))*I(:,1);
rho_fock = psi_fock*psi_fock';
g2_fock = trace(a'*a'*a*a*rho_fock)/trace(a'*a*rho_fock)^2;
g2_fock_check = (p_number - 1)/p_number;

%% 
rho_s = p1*rho_fock + p2*rho_coh + p3*rho_th;


g2_s_an = (p1*p_number*(p_number-1) + p2*alpha^4 + 2*p3*nth^2)/(p1*p_number + p2*alpha^2 + p3*nth)^2;
g2_s_nu = trace(a'*a'*a*a*rho_s)/trace(a'*a*rho_s)^2;

end