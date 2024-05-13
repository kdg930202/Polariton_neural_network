clc;
clearvars;
close all


d = 4; %dimension of the annihilation and creation operator
J = 1;
W = rand(1,2);
X = 0.5;
K = 1;
P = 0.2; %Incohernet pumping
wc = 0;
wl = 10;
swtich_H_I = 1;

dt = 0.1;
T = 0:dt:500;

Inten = 20;
E = linspace(0,Inten, length(T));

a = diag(sqrt(1:d-1),1); %annihilation operator
I_a = eye(d);


%Thermal state
nth = 0.5;
Rho_th = 0;


for i=0:d-1
    Rho_th = Rho_th+(nth^i/(1+nth)^(i+1))*I_a(:,i+1)*I_a(:,i+1)';
end




rho = Rho_th;
%%

for t=1:length(T)


    full(t) = trace(rho*a);
    g2(t) = trace(rho*a'*a'*a*a)/trace(rho*a'*a);


    K1 = -1i*(wc - wl)*(a'*a*rho - rho*a'*a) - 1i*X*(a'*a'*a*a*rho - rho*a'*a'*a*a) ...
         + (E(t)*a' - E(t)*a)*rho - rho*(E(t)*a' - E(t)*a) ...
         + K*(2*a*rho*a' - rho*a'*a - a'*a*rho + 2*nth*(a*rho*a' - rho*a*a' - a'*a*rho + a'*rho*a));

    rho1 = rho + 0.5*dt*K1;

    K2 = -1i*(wc - wl)*(a'*a*rho1 - rho1*a'*a) - 1i*X*(a'*a'*a*a*rho1 - rho1*a'*a'*a*a) ...
         + (E(t)*a' - E(t)*a)*rho1 - rho1*(E(t)*a' - E(t)*a) ...
         + K*(2*a*rho1*a' - rho1*a'*a - a'*a*rho1 + 2*nth*(a*rho1*a' - rho1*a*a' - a'*a*rho1 + a'*rho1*a));

    rho2 = rho + 0.5*dt*K2;
    
    K3 = -1i*(wc - wl)*(a'*a*rho2 - rho2*a'*a) - 1i*X*(a'*a'*a*a*rho2 - rho2*a'*a'*a*a) ...
         + (E(t)*a' - E(t)*a)*rho2 - rho2*(E(t)*a' - E(t)*a) ...
         + K*(2*a*rho2*a' - rho2*a'*a - a'*a*rho2 + 2*nth*(a*rho2*a' - rho2*a*a' - a'*a*rho2 + a'*rho2*a));

    rho3 = rho + dt*K3;

    K4 = -1i*(wc - wl)*(a'*a*rho3 - rho3*a'*a) - 1i*X*(a'*a'*a*a*rho3 - rho3*a'*a'*a*a) ...
         + (E(t)*a' - E(t)*a)*rho3 - rho3*(E(t)*a' - E(t)*a) ...
         + K*(2*a*rho3*a' - rho3*a'*a - a'*a*rho3 + 2*nth*(a*rho3*a' - rho3*a*a' - a'*a*rho3 + a'*rho3*a));

    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end


figure()
plot(E,abs(full))
