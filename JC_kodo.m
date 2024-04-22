clearvars
close all

d=5; %dimension of the annihilation and creation operator
hbar = 1;
wa = 1; %Atomic frequency
wf = 1; %Field frequency
g = 0.1; %Interaction constant 0.1
k = 0.05; %decay rate of cavity
dt = 0.1;
T = 0:dt:150;

a = diag(sqrt(1:d-1),1); %annihilation operator
If = eye(d);
Ia = eye(2);
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];
gs = [0;1]; %Ground state
es = [1;0]; %Excited state

%coherent state
alpha = 0.6;
coh = 0;

for i=0:d-1
    coh = coh+exp(-(alpha^2)/2)*alpha^i/sqrt(prod(1:i))*If(:,i+1);
end

g_coh = coh'*(a'*a'*a*a)*coh/(coh'*(a'*a)*coh);

%Thermal state
nth = 0.2;
the = 0;
for i=0:d-1
    the = the+(nth/(1+nth))^i*If(:,i+1)*If(:,i+1)';
end


H_atom = 0.5*hbar * wa * kron(sig_z,If);
H_field = hbar * wf * kron(Ia,a'*a);
H_inter = hbar*g*(kron(sig_m,a')+kron(sig_p,a));
H = H_atom + H_field + H_inter;

Ini_atom = es;
Ini_cavi = coh;%If(:,1)

psi = kron(Ini_atom, Ini_cavi); %Initial state=>atom in excited state and cavity in vaccum
rho = psi*psi'; %Initial density matrix

A = kron(Ia, a);

for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;
    K1 = -1i*(H*rho - rho*H) + 0.5*k*(2*A*rho*A' - rho*A'*A - A'*A*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + 0.5*k*(2*A*rho1*A' - rho1*A'*A - A'*A*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + 0.5*k*(2*A*rho2*A' - rho2*A'*A - A'*A*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + 0.5*k*(2*A*rho3*A' - rho3*A'*A - A'*A*rho3);
    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end

% hold on
% plot(T,Pe)
% plot(T,abs(correlation))
