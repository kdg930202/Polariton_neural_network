clc;
clearvars;

e1 = 0.5;
e2 = 0.5;

gamma1 = 2;
gamma2 = 2;
%%
E = 2 / sqrt(e1 * gamma1);

d=2; %dimension of the annihilation and creation operator

dt = 0.1;
T = 0:dt:150;

a = diag(sqrt(1:d-1),1); %annihilation operator
% If = eye(d);
Iden = eye(2);
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];
gs = [0;1]; %Ground state
es = [1;0]; %Excited state




sig_m_1 = kron(Iden, sig_m);
sig_m_2 = kron(sig_m, Iden);

sig_p_1 = kron(Iden, sig_p);
sig_p_2 = kron(sig_p, Iden);


H = -1i*sqrt(e1*gamma1)*(E*sig_m_1 - conj(E)*sig_m_1');
L1 = gamma1*(kron(sig_m_1', sig_m_1) - 0.5*kron(Iden, sig_m_1'*sig_m_1) - 0.5*kron(sig_m_1'*sig_m_1, Iden));
L2 = gamma2*(kron(sig_m_2', sig_m_2) - 0.5*kron(Iden, sig_m_2'*sig_m_2) - 0.5*kron(sig_m_2'*sig_m_2, Iden));
L3 = sqrt((1-e1)*(1-e1)*gamma1*gamma2) * (kron(sig_p_2, sig_m_1) - kron(sig_p_1, sig_m_2));
L4 = sqrt((1-e1)*(1-e1)*gamma1*gamma2) * (kron(sig_p_2, sig_m_1) - kron(sig_p_1, sig_m_2));

for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    % correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;
    K1 = -1i*(H*rho - rho*H) + gamma1*(a*rho1);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + 0.5*k*(2*A*rho1*A' - rho1*A'*A - A'*A*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + 0.5*k*(2*A*rho2*A' - rho2*A'*A - A'*A*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + 0.5*k*(2*A*rho3*A' - rho3*A'*A - A'*A*rho3);
    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end
