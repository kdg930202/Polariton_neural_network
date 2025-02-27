function [n_s, g2_s] = g2_ns_gs(phi,theta)
% clearvars
% close all
% clc

n1_selected = 1;
n2_selected = 1;
J = 1;
dt = 0.1;
T = 0:dt:20;
T = round(T,1);
time = 10;
tau = 1.5;
gamma = 1;
P = 0.1; %Incohernet pumping


% theta = pi/4;
% phi = pi/4;
p1 = cos(theta)^2;
p2 = sin(theta)^2 * cos(phi)^2;
p3 = sin(theta)^2 * sin(phi)^2;


d = 6; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
p_number = 1;
I = eye(d);
alpha = 1;
% TD = 100000000000;
TD = 1;

W = 1*rand(1,2);
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];

gs = [0;1]; %Ground state
es = [1;0]; %Excited state
Ini = 1/sqrt(2)*[1;1];

I_a = eye(d);
I_b = eye(2);
b1 = kron(sig_m, kron(I_b, I_a));
% b1 = kron(I_a,kron(sig_m, I_b));
b2 = kron(I_b, kron(sig_m, I_a));
% b2 = kron(I_a, kron(I_b, sig_m));
A = kron(I_b, kron(I_b,a));
% A = kron(a,kron(I_b,I_b));

H_R = J*(b1'*b2 + b2'*b1);

psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
rho_R = psi*psi';

%%

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
g2_s = trace(a'*a'*a*a*rho_s)/trace(a'*a*rho_s)^2;
n_s = trace(a'*a*rho_s);

end