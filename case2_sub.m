clearvars
close all


% m = 0:1:5;
r = linspace(0.2,1,1000);
% r=0;
% result_nu = zeros(length(phi), length(theta));
n_all = cell(1,length(r));
% n2_selected = cell(1,length(r));
% g2_all = zeros(length(m),length(r));

for i = 1:length(r)
        disp(i)
        % [g2_m(i), n1_selected{i}, n2_selected{i}] = case2_data_gen(r(i));
        [n_all{i}, g2_m(i)] = case2_subt(2,r(i));
        
end

function [n_all, g2_m] = case2_subt(m,r)
% function g2_m = case2_displaced(m,r)


J = 1;
dt = 0.1;
T = 0:dt:20;
% t = 0:dt:T;
time = 10;
tau = 1.5;
gamma = 1;
P = 0.1; %Incohernet pumping
% r = 0.1;
% m = 5;
TD = 100000000;

d = 20; %dimension of the annihilation and creation operator
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
g2_m = trace(a'*a'*a*a*rho_m)/trace(a'*a*rho_m)^2;
n_m = trace(a'*a*rho_m);


rho = kron(rho_m, rho_R);

for t=1:length(T)
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




select_step = 0.5;
select_points = time:select_step:(time+4-0.5);

for i =1:length(select_points)
    n1_selected(i) = n1(find(T==select_points(i)));
    n2_selected(i) = n2(find(T==select_points(i)));
end

n_all = [n1_selected, n2_selected];



% 
% figure()
% subplot(1,2,1)
% plot(T,abs(n1))
% hold on
% scatter(select_points, n1_selected)
% subplot(1,2,2)
% plot(T,abs(n2))
% hold on
% scatter(select_points, n2_selected)
end