clc;
clearvars;
close all

%%%
d = 4; %dimension of the annihilation and creation operator
dt = 0.01;
T = 0:dt:10;
ome_a = 1;
ome_b = 1;
U = 1;
g = 1;
OME_a = 1;
OME_b = 1;
hbar = 1;
gamma = 1;

a = diag(sqrt(1:d-1),1); %annihilation operator
I_a = eye(d);
I_b = eye(2);
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];
gs = [0;1]; %Ground state
es = [1;0]; %Excited state

A = kron(a, I_b);
B = kron(I_a, sig_m);


H = hbar*ome_a/(2*pi)*A'*A + hbar*ome_a/(2*pi)*B'*B + U/2*B'*B'*B*B ...
    +hbar*g/(2*pi)*(A'*B + A*B');


%%
rho = zeros(length(A(:,end)),length(A(:,end)));
rho(1,1) = 1;


% tau = 100;
% time = 500;
for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    % correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;

    K1 = -1i*(H*rho - rho*H) + gamma/2*(2*A*rho*A' - rho*A'*A - A'*A*rho) ...
                             + gamma/2*(2*B*rho*B' - rho*B'*B - B'*B*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho - rho*H) + gamma/2*(2*A*rho1*A' - rho1*A'*A - A'*A*rho1) ...
                             + gamma/2*(2*B*rho1*B' - rho1*B'*B - B'*B*rho1);

    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho - rho*H) + gamma/2*(2*A*rho2*A' - rho2*A'*A - A'*A*rho2) ...
                             + gamma/2*(2*B*rho2*B' - rho2*B'*B - B'*B*rho2);

    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho - rho*H) + gamma/2*(2*A*rho3*A' - rho3*A'*A - A'*A*rho3) ...
                             + gamma/2*(2*B*rho3*B' - rho3*B'*B - B'*B*rho3);

    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
    n(t) = trace(rho*B'*B);
end


figure()
plot(T,abs(n))


trace(rho)
g2 = trace(A'*A'*A*A*rho)/trace(A'*A*rho)^2;