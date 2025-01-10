clearvars
close all
clc

d = 6; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
nth=0.5;
gamma = 1;
alpha = 1;
p_number = 1;
a1 = a;
w = 1;
% T = 100;
dt = 0.1;
T = 0:dt:20;

theta = pi/4;
phi = pi/4;
p1 = cos(theta)^2;
p2 = sin(theta)^2 * cos(phi)^2;
p3 = sin(theta)^2 * sin(phi)^2;

Ia1 = eye(length(a1));
a2 = a;
Ia2 = eye(length(a2));
a3 = a;
Ia3 = eye(length(a3));
A1 = kron(a1, kron(Ia2,Ia3));
A2 = kron(Ia1, kron(a2,Ia3));
A3 = kron(Ia1, kron(Ia2,a3));

rho = zeros(length(a3),length(a3));
rho(1,1)=1;


for t=1:length(T)

    n_a(t) = trace(rho*a3'*a3); %this should be constant when u add cascade coupling term
    
    % H = H_R;
    K1 = -1i*w*(a3'*a3*rho - rho*a3'*a3) + gamma/2*(nth+1)*(2*a3*rho*a3' - a3'*a3*rho - rho*a3'*a3) ...
                                         + gamma/2*nth*(2*a3'*rho*a3 - a3*a3'*rho - rho*a3*a3');

    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*w*(a3'*a3*rho1 - rho1*a3'*a3) + gamma/2*(nth+1)*(2*a3*rho1*a3' - a3'*a3*rho1 - rho1*a3'*a3) ...
                                           + gamma/2*nth*(2*a3'*rho1*a3 - a3*a3'*rho1 - rho1*a3*a3');

    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*w*(a3'*a3*rho2 - rho2*a3'*a3) + gamma/2*(nth+1)*(2*a3*rho2*a3' - a3'*a3*rho2 - rho2*a3'*a3) ...
                                           + gamma/2*nth*(2*a3'*rho2*a3 - a3*a3'*rho2 - rho2*a3*a3');

    rho3 = rho + dt*K3;
    K4 = -1i*w*(a3'*a3*rho3 - rho3*a3'*a3) + gamma/2*(nth+1)*(2*a3*rho3*a3' - a3'*a3*rho3 - rho3*a3'*a3) ...
                                           + gamma/2*nth*(2*a3'*rho3*a3 - a3*a3'*rho3 - rho3*a3*a3');


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end

% trace(rho)
% trace(A3'*A3*rho)
figure()
plot(T,abs(n_a))
% ylim([0.3,0.5])
rho3 = rho;

%%
vacu = zeros(1,length(a2))';
vacu(1)=1;
% rho22 = expm(alpha*A2' - conj(alpha)*A2)*vacu*expm(alpha*A2' - conj(alpha)*A2);
psi_coh = expm(alpha*a2' - conj(alpha)*a2)*vacu;
rho2 = psi_coh*psi_coh';
% [trace(rho2), trace(A2'*A2*rho2), trace(A2'*A2'*A2*A2*rho2)/trace(A2'*A2*rho2)^2]
% rho2 = rho2*vacu'

%% 

II = zeros(1,length(a1))';
II(1) = 1;
psi_fock = (a1')^p_number/sqrt(prod(1:p_number))*II(:,1);
rho1 = psi_fock*psi_fock';
% g2_fock = trace(a3'*a3'*a3*a3*rho_fock)/trace(a3'*a3*rho_fock)^2;
% g2_fock_check = (p_number - 1)/p_number;

%%
tra = [trace(rho1), trace(rho2), trace(rho3)];
particle_n = [trace(a1'*a1*rho1), trace(a2'*a2*rho2), trace(a3'*a3*rho3)];
fprintf('trace_rho : %.2f %.2f %.2f \n',tra)
fprintf('particle_num: %.2f %.2f %.2f \n',particle_n)
fprintf('n: %.2f, alpha^2: %.2f, nth: %.2f \n',p_number, alpha^2, nth)


rho_s = kron(rho1,kron(rho2,rho3));
A = p1*A1 + p2*A2 + p3*A3;

g2 = trace(A'*A'*A*A*rho_s)/trace(A'*A*rho_s)^2;