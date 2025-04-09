clearvars
close all
clc


rr = 0:0.01:1;

for i = 1:length(rr)
    [g2_nu(i), g2_an(i), n_m(i), n_an(i)] = case2__(rr(i));
end

plot(rr,g2_nu,LineWidth=4)
hold on
plot(rr,g2_an,'y--',LineWidth=2)
text(0.8, 1, 'g2',FontSize=20)

% legend('Numerics','Analytics',fontsize=10)


plot(rr,n_m,LineWidth=4)
hold on
plot(rr,n_an,'y--',LineWidth=2)
text(0.5, 2.5, 'n',FontSize=20)
% legend('Numerics','Analytics',fontsize=10)

function [g2_nu, g2_an, n_m, n_an] = case2__(r)

%%
% clearvars
J = 1;
dt = 0.1;
T = 0:dt:20;
% t = 0:dt:T;
time = 10;
tau = 1.5;
gamma = 1;
P = 0.1; %Incohernet pumping
% r = 0.1;
m=1;
TD = 100000000;

d = 50; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
p_number = 1;
I = eye(d);
alpha = 1;


W = [0.5,0.5];
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];

gs = [0;1]; %Ground state
es = [1;0]; %Excited state
Ini = 1/sqrt(2)*[1;1];

I_a = eye(d);
I_b = eye(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sequence of the operators : source -> reservoir
% then the sequence of the initial density matrix 
% should be also rho_a(X)rho_b
A = kron(a, kron(I_b,I_b));
b1 = kron(I_a, kron(sig_m, I_b));
b2 = kron(I_a, kron(I_b, sig_m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_R = J*(b1'*b2 + b2'*b1);

psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
rho_R = psi*psi';


%%
S = expm(0.5*(r'*a*a - r*a'*a')); 
% D = expm(r*a'-r'*a);
vacc = zeros(length(a(:,1)),length(a(:,1)));
vacc(1,1) = 1;

adm = I_a;
% m = 0;
for i=1:m
    adm = a'*adm;
end
% psi_m = adm*D*S*vacc;
psi_m = adm*S*vacc;
psi_m = psi_m/norm(psi_m);

rho_m = psi_m*psi_m';
g2_nu = trace(a'*a'*a*a*rho_m)/trace(a'*a*rho_m)^2;
n_m = trace(a'*a*rho_m);
sinh_r = sinh(r);
% g2_an = (sinh_r^4 + 2 * sinh_r^2) / (sinh_r^2 + 1)^2;
n_m = trace(a'*a*rho_m);
n_an = 3*cosh(r)^2 - 2;
g2_an = 3*(3+2*tanh(r)^2)/(coth(r)-tanh(r))^2/n_an^2;
% [n_m, sinh(r)^2, trace(rho_m)]


end