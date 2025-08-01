%Hamiltonian:
%H=(Delta1)*(a1^t*a1)+(Delta2)*(a2^t*a2)+
%hbar*J1*(a1^t*a2+a2^t*a1)+
%0.5*U1*(a1^t*a1^t*a1*a1)+
%hbar*(E1)*(a1^t+a1)+hbar*(E2)*(a2^t+a2)
clc;
clearvars;
close all;
%----constants
gamma1=1;% all the parametrs are scaled in gamma1
gamma2=0.1;
% gamma2=0.05;
% dDel = 0.002;
dDel = 0.5;
% Delta1_vector = -1:dDel:1-dDel;
Delta1_vector = 0;
%%
Delta2 = 1.2;
% Delta2=10;
% E1_vector=0.01:0.02:0.5;% keep it in the low pumping regime
% E1 = omega_a, E2 = omega_b
E1_vector = 0.25;
phi=pi/2;
E2=0.4*exp(1i*phi);
J=1.6;%in the strong coupling regime it should be larger than decay rate
U=1;%
%------
Na1 = 6;
Na2 = 6;
dima1 = Na1 + 1;
dima2 = Na2 + 1;
dimtot=dima1*dima2;
Ia1 = speye(dima1); %identity on mode a subspace
Ia2 = speye(dima2); %identity on cascade system subspace
Itot = speye(dimtot); %identity on the complete Hilbert space
% "a" mode matric and "b" modde matrix
a1 = spdiags(sqrt(0:Na1)',1,dima1,dima1);
a2=spdiags(sqrt(0:Na2)',1,dima2,dima2);
a1 = kron(a1,Ia2);
a2 = kron(Ia1,a2);
%-----main
rho_s = zeros(dimtot, dimtot, length(Delta1_vector), length(E1_vector));
for k=1:length(E1_vector)
    E1=E1_vector(k);
    for j=1:length(Delta1_vector)
        Delta1=Delta1_vector(j);
        H1=(Delta1)*(a1'*a1)+(Delta2)*(a2'*a2);
        H2=J*(a1'*a2+a2'*a1);
        H3=conj(E1)*a1'+E1*a1+conj(E2)*a2'+E2*a2;
        H4=U*(a1'*a1'*a1*a1);
        H=H1+H2+H3+H4;
        La1=(1/2)*gamma1*(2*kron(conj(a1),a1)-kron(Itot,a1'*a1)-kron(a1.'*conj(a1),Itot));
        La2=(1/2)*gamma2*(2*kron(conj(a2),a2)-kron(Itot,a2'*a2)-kron(a2.'*conj(a2),Itot));
        L=-1i*kron(Itot,H)+1i*kron(H.',Itot)+La1+La2;
    %Firts method
     [rhoS1,lambda0] = eigs(L,1,'sm'); 
     eigen0 = lambda0;
     rhoS1 = reshape(rhoS1,dimtot,dimtot);
     rhoS1 = rhoS1/trace(rhoS1);
     rho_s(:, :, j,k)=rhoS1;
     %pop(j)= trace(a1'*a1*rhoS1);
     g2_a(k,j) = trace((a1'*a1*a1'*a1-a1'*a1)*rhoS1)/(trace(a1'*a1*rhoS1))^2;
     n_a(k,j) = trace(a1'*a1*rhoS1);
     % g2_b(k,j)=trace((a2'*a2*a2'*a2-a2'*a2)*rhoS1)/(trace(a2'*a2*rhoS1))^2;
     disp(j)
    end
end

%%
plot(Delta1_vector, abs(g2_a))
hold on
plot(Delta1_vector, abs(n_a))
xlabel('Delta_1')
ylabel('g2')
% legend(['g2','population of the source'])
%%
W = [0.5,0.5];
ns = cell(1,length(Delta1_vector));

for i=1:length(Delta1_vector)
    % disp('reservoir : ',i)
    fprintf('reservoir step %d\n',i);
    ns{i} = case3_(W,rho_s(:,:,i));
end


function ns = case3_(W,rho_s)
dt = 0.1;
T = 0:dt:10-dt;
T = round(T,1);
time = 6;
tau = 1;
J = 1;
trun_de = 0;
% W = 1*rand(1,2);
pump = 3;
gamma = 1;
P = 0.1;
TD = 100000000;



d = 7; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
a1 = a;
a2 = a;
sig_z = [1,0;0,-1];
sig_m = [0,0;1,0];
sig_p = [0,1;0,0];

gs = [0;1]; %Ground state
es = [1;0]; %Excited state
% Ini = 1/sqrt(2)*[1;1];
% 
% A = kron(a, kron(I_b,I_b));
% b1 = kron(I_a, kron(sig_m, I_b));
% b2 = kron(I_a, kron(I_b, sig_m));

Ia1 = eye(d); %identity on mode a subspace
Ia2 = eye(d); %identity on cascade system subspace
I_b = eye(2);

A1 = sparse(kron(kron(a1, Ia1), kron(I_b,I_b)));
% A2 = kron(kron(Ia1, a2), kron(I_b,I_b));
b1 = sparse(kron(kron(Ia1, Ia1), kron(sig_m,I_b)));
b2 = sparse(kron(kron(Ia1, Ia2), kron(I_b,sig_m)));

H_R = J*(b1'*b2 + b2'*b1);


psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
rho_R = sparse(psi*psi');


rho = kron(sparse(rho_s), rho_R);

for t=1:length(T)
    % disp(t)
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    % n_a(t) = trace(rho*A1'*A1); %this should be constant when u add cascade coupling term

    H = H_R;
    K1 = -1i*(H*rho - rho*H) + pump*(W(1)*(A1*rho*b1' - b1'*A1*rho + b1*rho*A1' - rho*A1'*b1) ...
                             + W(2)*(A1*rho*b2' - b2'*A1*rho + b2*rho*A1' - rho*A1'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                             + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho)...
                             + trun_de*1/(2*TD)*(2*b1'*b1*rho*b1'*b1 - rho*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho ...
                                       + 2*b2'*b2*rho*b2'*b2 - rho*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + pump*(W(1)*(A1*rho1*b1' - b1'*A1*rho1 + b1*rho1*A1' - rho1*A1'*b1) ...
                               + W(2)*(A1*rho1*b2' - b2'*A1*rho1 + b2*rho1*A1' - rho1*A1'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1)...
                               + trun_de*1/(2*TD)*(2*b1'*b1*rho1*b1'*b1 - rho1*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho1 ...
                                         + 2*b2'*b2*rho1*b2'*b2 - rho1*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + pump*(W(1)*(A1*rho2*b1' - b1'*A1*rho2 + b1*rho2*A1' - rho2*A1'*b1) ...
                               + W(2)*(A1*rho2*b2' - b2'*A1*rho2 + b2*rho2*A1' - rho2*A1'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                               + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                               + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                               + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2)...
                               + trun_de*1/(2*TD)*(2*b1'*b1*rho2*b1'*b1 - rho2*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho2 ...
                                         + 2*b2'*b2*rho2*b2'*b2 - rho2*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + pump*(W(1)*(A1*rho3*b1' - b1'*A1*rho3 + b1*rho3*A1' - rho3*A1'*b1) ...
                               + W(2)*(A1*rho3*b2' - b2'*A1*rho3 + b2*rho3*A1' - rho3*A1'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3)...
                               + trun_de*1/(2*TD)*(2*b1'*b1*rho3*b1'*b1 - rho3*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho3 ...
                                         + 2*b2'*b2*rho3*b2'*b2 - rho3*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho3);;


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end


select_step = dt;
sel_ini = 6+select_step;
sel_fin = 8;

select_points = sel_ini:select_step:sel_fin;
select_points = round(select_points,1);
for i =1:length(select_points)
    n1_selected(i) = n1(find(T==select_points(i)));
    n2_selected(i) = n2(find(T==select_points(i)));
end

ns = [abs(n1_selected), abs(n2_selected)];


figure(2)
plot(T,abs(n1))
hold on
plot([6,10],[0.0907,0.0907],'k--')
xlabel('Time, ${t\gamma/\hbar}$','interpreter','latex', 'FontSize',20)
ylabel('Reservoir density','FontSize',20)
scatter(select_points,abs(n1_selected),'r')


end