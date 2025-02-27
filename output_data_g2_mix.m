clc;
clear;
close all;
%-------------------RF
Delta_s=-2;
gamma_s=1;
Omega_s=0.8;
%------------------Laser
dDel = 0.005;
Delta_a_vector=-2.5:dDel:2.5-dDel;
gamma_a=1;
Omega_a=0.3;
%------------
theta=pi/4;



for i=1:length(Delta_a_vector)
    Delta_a=Delta_a_vector(i);
    [O2,g2_O1,g2_O2,rho_out,alpha_num,alpha_ana,pop_s,pop_a,n_O1,n_O2]=g2_mix(Delta_s,gamma_s,Omega_s,Delta_a,gamma_a,Omega_a,theta);
    g2_O1_4saving(i)=g2_O1;
    g2_O2_4saving(i)=g2_O2;
    rho_out_4saving(:,:,i)=rho_out;

    %density in laser |alpha|^2
    pop_laser_4saving(i)=pop_a;
    %density in RF
    pop_RF_4saving(i)=pop_s;
    %density of O1 out put
    pop_O1_4saving(i)=n_O1;
    %density of O2 out put
    pop_O2_4saving(i)=n_O2;

    %to check if N is selected correctly
    test_N(i)=abs(alpha_num)^2-abs(alpha_ana)^2;
end

% semilogy(Delta_a_vector,real(g2_O1_4saving),'linewidth',2);
% hold on
figure
semilogy(Delta_a_vector,real(g2_O2_4saving),'linewidth',2);
% xlabel('\Delat_a/ \gamma_a')
ylabel('g^{(2)}_{O_2}')
% legend('g'g^{(2)}_{O_2}')

figure;
 % plot(Delta_a_vector,real(pop_O1_4saving),'linewidth',2);
hold on
plot(Delta_a_vector,real(pop_O2_4saving),'linewidth',2);
% xlabel('\Delat_a/\gamma_a')
ylabel('n_{O_2}')
% legend('n_{O_1}','n_{O_2}')

% figure;
% we select N to make this almost zero
% plot(Delta_a_vector,test_N)
% hold on;
% plot(Delta_a_vector,pop_O1_4saving+pop_O2_4saving)

%%

gs = [0;1]; %Ground state
psi = kron(gs,gs); %reservoir 1 : excited, reservoir 2 : ground
rho_R = sparse(psi*psi');

dt = 0.1;
T = 0:dt:10-dt;
T = round(T,1);

select_step = dt;
sel_ini = 6+select_step;
sel_fin = 8;

select_points = sel_ini:select_step:sel_fin;
select_points = round(select_points,1);

all_n1 = NaN(length(Delta_a_vector), length(T));
all_ns = NaN(length(Delta_a_vector), 2*length(select_points));

for i=1:length(Delta_a_vector)
    disp(i)
    which_rho = rho_out_4saving(:,:,i);
    rho = kron(sparse(which_rho), rho_R);
    [n1, ns] = case4_reservoir(dt,T,select_points,O2,rho);
    all_n1(i,:) = abs(n1);
    all_ns(i,:) = abs(ns);
end


% for i = 1:length(Delta_a_vector)
%     figure(3)
%     plot(T,all_n1(i,:))
%     hold on
%     scatter(select_points,all_ns(i,1:length(select_points)))
% end

function [n1, ns] = case4_reservoir(dt,T,select_points,O2,rho)
time = 6;
tau = 1;
J = 1;
trun_de = 0;
% W = 1*rand(1,2);
W = [0.5,0.1];
pump = 3;
gamma = 1;
P = 0.1;
TD = 1;

sig_m = [0,0;1,0];

I_b = eye(2);
IO2 = eye(length(O2));
O2 = sparse(kron(O2, kron(I_b,I_b)));
% A2 = kron(kron(IO2, a2), kron(I_b,I_b));
b1 = sparse(kron(IO2, kron(sig_m,I_b)));
b2 = sparse(kron(IO2, kron(I_b,sig_m)));

H_R = J*(b1'*b2 + b2'*b1);

for t=1:length(T)
    % disp(t)
    n1(t) = trace(rho*b1'*b1);
    n2(t) = trace(rho*b2'*b2);
    % n_a(t) = trace(rho*O2'*O2); %this should be constant when u add cascade coupling term

    H = H_R;
    K1 = -1i*(H*rho - rho*H) + pump*(W(1)*(O2*rho*b1' - b1'*O2*rho + b1*rho*O2' - rho*O2'*b1) ...
                             + W(2)*(O2*rho*b2' - b2'*O2*rho + b2*rho*O2' - rho*O2'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                             + gamma/2*(2*b1*rho*b1' - rho*b1'*b1 - b1'*b1*rho) ...
                             + gamma/2*(2*b2*rho*b2' - rho*b2'*b2 - b2'*b2*rho) ...
                             + P/2*(2*b1'*rho*b1 - rho*b1*b1' - b1*b1'*rho)...
                             + P/2*(2*b2'*rho*b2 - rho*b2*b2' - b2*b2'*rho);
                             % + trun_de*1/(2*TD)*(2*b1'*b1*rho*b1'*b1 - rho*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho ...
                             %           + 2*b2'*b2*rho*b2'*b2 - rho*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho);
    rho1 = rho + 0.5*dt*K1;
    K2 = -1i*(H*rho1 - rho1*H) + pump*(W(1)*(O2*rho1*b1' - b1'*O2*rho1 + b1*rho1*O2' - rho1*O2'*b1) ...
                               + W(2)*(O2*rho1*b2' - b2'*O2*rho1 + b2*rho1*O2' - rho1*O2'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho1*b1' - rho1*b1'*b1 - b1'*b1*rho1) ...
                               + gamma/2*(2*b2*rho1*b2' - rho1*b2'*b2 - b2'*b2*rho1) ...
                               + P/2*(2*b1'*rho1*b1 - rho1*b1*b1' - b1*b1'*rho1) ...
                               + P/2*(2*b2'*rho1*b2 - rho1*b2*b2' - b2*b2'*rho1);
                               % + trun_de*1/(2*TD)*(2*b1'*b1*rho1*b1'*b1 - rho1*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho1 ...
                                         % + 2*b2'*b2*rho1*b2'*b2 - rho1*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho1);
    rho2 = rho + 0.5*dt*K2;
    K3 = -1i*(H*rho2 - rho2*H) + pump*(W(1)*(O2*rho2*b1' - b1'*O2*rho2 + b1*rho2*O2' - rho2*O2'*b1) ...
                               + W(2)*(O2*rho2*b2' - b2'*O2*rho2 + b2*rho2*O2' - rho2*O2'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho2*b1' - rho2*b1'*b1 - b1'*b1*rho2) ...
                               + gamma/2*(2*b2*rho2*b2' - rho2*b2'*b2 - b2'*b2*rho2) ...
                               + P/2*(2*b1'*rho2*b1 - rho2*b1*b1' - b1*b1'*rho2) ...
                               + P/2*(2*b2'*rho2*b2 - rho2*b2*b2' - b2*b2'*rho2);
                               % + trun_de*1/(2*TD)*(2*b1'*b1*rho2*b1'*b1 - rho2*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho2 ...
                                         % + 2*b2'*b2*rho2*b2'*b2 - rho2*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho2);
    rho3 = rho + dt*K3;
    K4 = -1i*(H*rho3 - rho3*H) + pump*(W(1)*(O2*rho3*b1' - b1'*O2*rho3 + b1*rho3*O2' - rho3*O2'*b1) ...
                               + W(2)*(O2*rho3*b2' - b2'*O2*rho3 + b2*rho3*O2' - rho3*O2'*b2))*(T(t)>=time).*(T(t)<time+tau) ...
                               + gamma/2*(2*b1*rho3*b1' - rho3*b1'*b1 - b1'*b1*rho3) ...
                               + gamma/2*(2*b2*rho3*b2' - rho3*b2'*b2 - b2'*b2*rho3) ...
                               + P/2*(2*b1'*rho3*b1 - rho3*b1*b1' - b1*b1'*rho3) ...
                               + P/2*(2*b2'*rho3*b2 - rho3*b2*b2' - b2*b2'*rho3);
                               % + trun_de*1/(2*TD)*(2*b1'*b1*rho3*b1'*b1 - rho3*b1'*b1*b1'*b1 - b1'*b1*b1'*b1*rho3 ...
                               %           + 2*b2'*b2*rho3*b2'*b2 - rho3*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho3);


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);
end



for j = 1:length(select_points)
    n1_selected(j) = n1(find(T==select_points(j)));
    n2_selected(j) = n2(find(T==select_points(j)));
end

ns = [abs(n1_selected), abs(n2_selected)];
end
    % for j = 1:length(select_points)
    %     n1_selected(j) = n1(find(T==select_points(j)));
    %     n2_selected(j) = n2(find(T==select_points(j)));
    % end
    % big_n1(i,:) = abs(n1);
    % ns(i,:) = [abs(n1_selected),abs(n2_selected)];



%%
% figure()
% plot(T,abs(n1))
% hold on
% plot(T,abs(n2))
% scatter(select_points, abs(n1_selected))
% scatter(select_points, abs(n2_selected))



