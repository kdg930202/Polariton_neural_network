%Time evolution of a two level system under coherent pumping
clc;
clear;
close all;
%-------------------Constants
Gamma=0.1;
Omega=0;
E1=0;
E2=1;

%-----------------Operators
sigma=[0 0;1 0];
sigma_z = sigma'*sigma - sigma*sigma';
I_sigma=eye(2);
%------------------Hamiltonian
% H=Omega*sigma'*sigma+E*(sigma'+sigma);
H = 0.5*(E2-E1)*sigma_z + 0.5*(E1+E2)*I_sigma + Omega*(sigma'+sigma);

%%
%-------------------Time
T_max=20;nkk=2;
dt=0.025/nkk;Time=0:dt:T_max;Nt=length(Time);
%-----------------Initial conditions
p1=zeros(2,1);p1(1)=1;
rho=p1*p1';
%----------------------Dynamics
n_s=zeros(1,Nt);
for i=1:Nt
    n_s1(i)=trace(rho*sigma'*sigma); %Population of excited state
    n_s2(i)=trace(rho*sigma*sigma'); %Population of ground state
    % n_d(i)=trace(rho*sigma_z);
    rho1=rho-1i*(H*rho-rho*H)*dt+...
        (0.5*Gamma*(sigma*rho*sigma'-sigma'*sigma*rho+sigma*rho*sigma'-rho*sigma'*sigma))*dt;
    rho2=0.5*(rho+rho1);
    rho=rho-1i*(H*rho2-rho2*H)*dt+...
        (0.5*Gamma*(sigma*rho2*sigma'-sigma'*sigma*rho2+sigma*rho2*sigma'-rho2*sigma'*sigma))*dt;
end
%-------------------Ploting the result
% figure()
% subplot(1,2,1)
% plot(Time,real(n_s1),'LineWidth',1.5)
% xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
% ylabel('$<\sigma^\dagger \sigma>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
% ax = gca;
% ax.FontSize = 20;
% pbaspect([1 1 1])
% subplot(1,2,2)
% plot(Time,real(n_s2),'LineWidth',1.5)
% xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
% ylabel('$<\sigma_z>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
% ax = gca;
% ax.FontSize = 20;
% pbaspect([1 1 1])

figure()
plot(Time, abs(n_s1))
hold on
plot(Time, abs(n_s2))
% plot(Time, real(n_s1) + real(n_s2))