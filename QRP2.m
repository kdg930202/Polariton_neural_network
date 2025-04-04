clc;
clearvars;
close all

%%%
d = 4; %dimension of the annihilation and creation operator
J = 1;
W = rand(1,2);
gamma = 1;
P = 0.2; %Incohernet pumping
swtich_H_I = 1;

dt = 0.1;
T = 0:dt:150;

a = diag(sqrt(1:d-1),1); %annihilation operator
I_a = eye(d);
I_b = eye(2);
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];
gs = [0;1]; %Ground state
es = [1;0]; %Excited state


b1 = kron(sig_m, kron(I_b, I_a));
b2 = kron(I_b, kron(sig_m, I_a));
A = kron(I_b, kron(I_b,a));
% a2 = kron(I_b, a);

H_R = J*(b1'*b2 + b2'*b1);
% H_I = W(1)*(A'*b1 + b1'*A) + W(2)*(A'*b2 + b2'*A );


% b11 = kron
psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
Rho_b = psi*psi';



%psi = [0;0;0;0.1];


%Thermal state
nth = 0.1;
Rho_th = 0;

for i=0:d-1
    Rho_th = Rho_th+(nth^i/(1+nth)^(i+1))*I_a(:,i+1)*I_a(:,i+1)';
    % Rho_th = Rho_th + (nth/(1+nth))^i * I_a(:,i+1)*I_a(:,i+1)';
end

for i=0:d-1
    Pn_th(i+1) = I_a(:,i+1)'*Rho_th*I_a(:,i+1);
end

% figure()
% bar(0:d-1, Pn_th, 'r')

% rho = kron(Rho_th, Rho_b);
rho = kron(Rho_b, Rho_th);
trace(rho*b1'*b1)
% trace(Rho_th) %should be 1
% trace(Rho_th*a'*a) %should be n_th
% trace(rho) %should be 1
% trace()


tau = 20;
time = 500;
for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    % correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    n_a(t) = trace(rho*A'*A); %this should be constant when u add cascade coupling term
    
    % H_I = W(1)*(A*rho*b1' - b1'*A*rho + b1*rho*A' - rho*A'*b1) ...
    %     + W(2)*(A*rho*b2' - b2'*A*rho + b2*rho*A' - rho*A'*b2);
    H = H_R;
    K1 = -1i*(H*rho - rho*H) + (W(1)*(A*rho*b1' - b1'*A*rho + b1*rho*A' - rho*A'*b1) ...
                             + W(2)*(A*rho*b2' - b2'*A*rho + b2*rho*A' - rho*A'*b2))*(t>=time).*(t<time+tau) ...
                             + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + (W(1)*(A*rho1*b1' - b1'*A*rho1 + b1*rho1*A' - rho1*A'*b1) ...
                               + W(2)*(A*rho1*b2' - b2'*A*rho1 + b2*rho1*A' - rho1*A'*b2))*(t>=time).*(t<time+tau) ...
                               + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + (W(1)*(A*rho2*b1' - b1'*A*rho2 + b1*rho2*A' - rho2*A'*b1) ...
                               + W(2)*(A*rho2*b2' - b2'*A*rho2 + b2*rho2*A' - rho2*A'*b2))*(t>=time).*(t<time+tau) ...
                               + gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                               + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                               + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                               + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + (W(1)*(A*rho3*b1' - b1'*A*rho3 + b1*rho3*A' - rho3*A'*b1) ...
                               + W(2)*(A*rho3*b2' - b2'*A*rho3 + b2*rho3*A' - rho3*A'*b2))*(t>=time).*(t<time+tau) ...
                               + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3);

    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end


figure()
plot(T,abs(n1))
hold on
% plot(T,abs(n2))
% 
figure()
plot(T,abs(n_a))