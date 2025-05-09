function [n_s, g2_s, n1_selected, n2_selected] = g2_EQ13(W, phi, theta)
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

% W = 1*rand(1,2);
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

rho = kron(rho_R, rho_s);

for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    % correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    n_a(t) = trace(rho*A'*A); %this should be constant when u add cascade coupling term

    H = H_R;
    K1 = -1i*(H*rho - rho*H) + (W(1)*(A*rho*b1' - b1'*A*rho + b1*rho*A' - rho*A'*b1) ...
                             + W(2)*(A*rho*b2' - b2'*A*rho + b2*rho*A' - rho*A'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                             + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho)...
                             + 1/(2*TD)*(2*b1'*b1*rho*b1'*b1 - rho*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho ...
                                       + 2*b2'*b2*rho*b2'*b2 - rho*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + (W(1)*(A*rho1*b1' - b1'*A*rho1 + b1*rho1*A' - rho1*A'*b1) ...
                               + W(2)*(A*rho1*b2' - b2'*A*rho1 + b2*rho1*A' - rho1*A'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1)...
                               + 1/(2*TD)*(2*b1'*b1*rho1*b1'*b1 - rho1*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho1 ...
                                         + 2*b2'*b2*rho1*b2'*b2 - rho1*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + (W(1)*(A*rho2*b1' - b1'*A*rho2 + b1*rho2*A' - rho2*A'*b1) ...
                               + W(2)*(A*rho2*b2' - b2'*A*rho2 + b2*rho2*A' - rho2*A'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                               + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                               + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                               + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2)...
                               + 1/(2*TD)*(2*b1'*b1*rho2*b1'*b1 - rho2*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho2 ...
                                         + 2*b2'*b2*rho2*b2'*b2 - rho2*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + (W(1)*(A*rho3*b1' - b1'*A*rho3 + b1*rho3*A' - rho3*A'*b1) ...
                               + W(2)*(A*rho3*b2' - b2'*A*rho3 + b2*rho3*A' - rho3*A'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3)...
                               + 1/(2*TD)*(2*b1'*b1*rho3*b1'*b1 - rho3*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho3 ...
                                         + 2*b2'*b2*rho3*b2'*b2 - rho3*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho3);;


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end

% figure()
% subplot(1,2,1)
% plot(T,abs(n1))
% subplot(1,2,2)
% plot(T,abs(n2))
% figure()
% plot(T,abs(n_a))
%%
select_step = 0.25;
select_points = time+select_step:select_step:(time+4-select_step);
select_points = round(select_points,1);
for i =1:length(select_points)
    n1_selected(i) = n1(find(T==select_points(i)));
    n2_selected(i) = n2(find(T==select_points(i)));
end


%%
% figure()
% plot(T, abs(n1))
% hold on
% scatter(select_points, abs(n1_selected))
end