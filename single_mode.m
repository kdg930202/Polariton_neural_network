clearvars

d = 6; %dimension of the annihilation and creation operator
J = 1;
aa = diag(sqrt(1:d-1),1); %annihilation operator
W = rand(1,2);
dt = 0.1;
T = 0:dt:150;
gamma = 1;
P = 0.2; %Incohernet pumping
swtich_H_I = 1;

gs = [0;1]; %Ground state
es = [1;0]; %Excited state

sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];

I_a = eye(d);
I_b = eye(2);
b1 = kron(sig_m, kron(I_b, I_a));
b2 = kron(I_b, kron(sig_m, I_a));
A = kron(I_b, kron(I_b,aa));

H_R = J*(b1'*b2 + b2'*b1);
H_I = W(1)*(A'*b1 + b1'*A) + W(2)*(A'*b2 + b2'*A );

% I_a = eye(d);
a1 = aa;

theta = 1.9605;
s = 0.1651;
phi = 0.5+pi/10;

alpha = s*sin(phi);
a = abs(alpha)*exp(1i*theta);

% Matrix exponential : expm
S = expm(0.5*(a'*a1*a1 - a*a1'*a1')); 


psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
Rho_b = psi*psi';


n_bar = s^2 * (cos(phi))^2 ;
% n_bar = 0.5;

for i=1:d
    rho_th(i) = (1/(1 + n_bar)) * (n_bar/(1 + n_bar))^(i-1);
    % Rho_th = Rho_th + rho_th * kron(I_a(:,i+1),I_a(:,i+1));
end

rho_th = diag(rho_th);
rho_sq_th = S*rho_th*conj(S);
rho = kron(Rho_b, rho_sq_th);
[trace(rho), trace(rho_sq_th),trace(rho*b1'*b1)]

% rho_th = diag(rho_th);
% [trace(rho_th),trace(rho_th*a1'*a1)] %should be 1 and n_bar
% 
% rho_sq_th = S*rho_th*S';
% % rho_sq_th = S*rho_th*conj(S);
% [trace(rho_sq_th),trace(rho_sq_th*a1'*a1)] %should be 1 and n_bar


tau = 10;
for t=1:length(T)
    % Pe(t) = kron(es, If(:,1))'*rho*kron(es, If(:,1));
    % correlation(t) = trace(rho*A'*A'*A*A)/trace(rho*A'*A)^2;
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    n_a(t) = trace(rho*A'*A); %this should be constant when u add cascade coupling term

    H = H_R + swtich_H_I * H_I.*(t>=500).*(t<500+tau);
    K1 = -1i*(H*rho - rho*H) + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H)+ gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                              + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                              + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                              + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3);

    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end


figure()
plot(T,abs(n1))
hold on
% plot(T,abs(n2))

