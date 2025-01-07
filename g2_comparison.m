function [g2_s_an, g2_s_nu] = g2_comparison(alpha, nth)

p1 = 1;
p2 = 1;
p3 = 1;
denom = sqrt(p1^2 + p2^2 + p3^2);
p1 = p1/denom;
p2 = p2/denom;
p3 = p3/denom;
% ps_sum = p1^2 + p2^2 + p3^2;
p1 = p1^2;
p2 = p2^2;
p3 = p3^2;
% ps_new_sum = p1 + p2 + p3; 
d = 12; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
p_number = 1;
I = eye(d);
% alpha = 0.1;

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
% nth = 0.5;
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