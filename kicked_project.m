%Time evolution of a two level system under coherent pumping
clc;
clear;
close all;
%-------------------Constants%
Delta = 1;
Omega = 1;
gamma = 1;
gamma_a = 0.5;
J = 1;
pump = 1;
P = 0.1;
%-----------------Operators
N = 6;
a = diag(sqrt(1:N-1),1);
Ia = eye(N);
%------------------Hamiltonian
% I_b = eye(2);
% sig_m = [0,0;1,0];
% W = [0.5,0.5];


% A2 = kron(kron(IA, a2), kron(I_b,I_b));


H=Delta*a'*a+Omega*(a'+a);





%-------------------Time
T_max = 20;
dt = 0.1;
Time = 0:dt:T_max-dt;
Time = round(Time,1);
which_T = 5;
tau = 1;
Nt = length(Time);
%-----------------Initial conditions
ket_gen = zeros(N,1); 
vec = ket_gen;
vec(1) = 1;
rho = vec*vec';
rho_s = cell(1,length(Time));




%----------------------Dynamics
n_s=zeros(1,Nt);
for i=1:Nt
    %second order
    rho_s{i} = rho;
    n_s(i)=trace(a'*a*rho);
    k1=-1i*(H*rho-rho*H)*dt+...
        dt*(0.5*gamma_a*(a*rho*a'-a'*a*rho+a*rho*a'-rho*a'*a));
    rho1=rho+0.5*k1;
    k2=-1i*(H*rho1-rho1*H)*dt+...
        dt*(0.5*gamma_a*(a*rho1*a'-a'*a*rho1+a*rho1*a'-rho1*a'*a));
    rho=rho+k2;

end
%-------------------Ploting the result
%%
figure;
% plot(Time,real(n_s),'LineWidth',1.5)
plot(abs(n_s),'LineWidth',1.5)
xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
ylabel('$<a^\dagger a>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
ax = gca;
ax.FontSize = 20;
pbaspect([1 1 1])

%%

ns_full = cell(1,length(Time));

for i=1:length(Time)
    disp(i)
    ns_full{i} = RC(rho_s{i}, Time, tau, which_T,a,Delta,Omega,J,pump, gamma,0, P, dt);
end


% 
% [n1,n2,select_points,ns] = RC(rho_s{25}, Time, tau, which_T,a,Delta,Omega,J,pump, gamma,gamma_a, P, dt);
% figure()
% plot(Time,abs(n1))
% hold on
% scatter(select_points,ns(1:length(ns)/2))


% function [n1, n2, select_points, ns] = RC(rho_s,Time,tau,which_T,a,Delta,Omega,J,pump, gamma,gamma_a, P, dt)
function ns = RC(rho_s,Time,tau,which_T,a,Delta,Omega,J,pump, gamma, gamma_a, P, dt)
Ia = eye(length(a));
I_b = eye(2);
sig_m = [0,0;1,0];
W = [0.5,0.5];
A = sparse(kron(a, kron(I_b,I_b)));
b1 = sparse(kron(Ia, kron(sig_m,I_b)));
b2 = sparse(kron(Ia, kron(I_b,sig_m)));
H_s = Delta*A'*A+Omega*(A'+A);
H_R = J*(b1'*b2 + b2'*b1);

gs = [0;1]; %Ground state
psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
rho_R = sparse(psi*psi');
rho = kron(rho_s,rho_R);

for t=1:length(Time)
    % disp(Time)
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    % n_a(Time) = trace(rho*A'*A); %this should be constant when u add cascade coupling term

    H = H_R+H_s;
    K1 = -1i*(H*rho - rho*H) + pump*(W(1)*(A*rho*b1' - b1'*A*rho + b1*rho*A' - rho*A'*b1) ...
                             + W(2)*(A*rho*b2' - b2'*A*rho + b2*rho*A' - rho*A'*b2))*(Time(t)>=which_T).*(Time(t)<which_T+tau) ...
                             + gamma_a/2*(2*A*rho*A' - A'*A*rho - rho*A'*A) ...
                             + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + pump*(W(1)*(A*rho1*b1' - b1'*A*rho1 + b1*rho1*A' - rho1*A'*b1) ...
                               + W(2)*(A*rho1*b2' - b2'*A*rho1 + b2*rho1*A' - rho1*A'*b2))*(Time(t)>=which_T).*(Time(t)<which_T+tau) ...
                               + gamma_a/2*(2*A*rho1*A' - A'*A*rho1 - rho1*A'*A) ...
                               + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + pump*(W(1)*(A*rho2*b1' - b1'*A*rho2 + b1*rho2*A' - rho2*A'*b1) ...
                               + W(2)*(A*rho2*b2' - b2'*A*rho2 + b2*rho2*A' - rho2*A'*b2))*(Time(t)>=which_T).*(Time(t)<which_T+tau) ...
                               + gamma_a/2*(2*A*rho2*A' - A'*A*rho2 - rho2*A'*A) ...
                               + gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                               + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                               + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                               + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + pump*(W(1)*(A*rho3*b1' - b1'*A*rho3 + b1*rho3*A' - rho3*A'*b1) ...
                               + W(2)*(A*rho3*b2' - b2'*A*rho3 + b2*rho3*A' - rho3*A'*b2))*(Time(t)>=which_T).*(Time(t)<which_T+tau) ...
                               + gamma_a/2*(2*A*rho3*A' - A'*A*rho3 - rho3*A'*A) ...
                               + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3);


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);

end

select_step = dt*4;
sel_ini = 5+select_step;
sel_fin = 8+select_step;

select_points = sel_ini:select_step:sel_fin;
select_points = round(select_points,1);
% display(select_points)
n1_selected = zeros(1,length(select_points));
n2_selected = zeros(1,length(select_points));
% display(n1(find(Time==select_points(1))))
for i =1:length(select_points)
    % display(n1(find(Time==select_points(i))))
    n1_selected(i) = n1(find(Time==select_points(i)));
    n2_selected(i) = n2(find(Time==select_points(i)));
end

ns = [abs(n1_selected), abs(n2_selected)];

end