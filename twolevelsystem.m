%Time evolution of a two level system under coherent pumping
clc;
clear;
close all;
%-------------------Constants
Gamma=5;
Omega=0;
E=1;
%-----------------Operators
sigma=[0 0;1 0];
I_sigma=eye(2);
%------------------Hamiltonian
H=Omega*sigma'*sigma+E*(sigma'+sigma);
%-------------------Time
T_max=11;nkk=2;
dt=0.025/nkk;Time=0:dt:T_max;Nt=length(Time);
%-----------------Initial conditions
p1=zeros(2,1);p1(2)=1;
rho=p1*p1';
%----------------------Dynamics
n_s=zeros(1,Nt);
for i=1:Nt
    n_s(i)=trace(rho*sigma'*sigma);
    rho1=rho-1i*(H*rho-rho*H)*dt+...
        (0.5*Gamma*(sigma*rho*sigma'-sigma'*sigma*rho+sigma*rho*sigma'-rho*sigma'*sigma))*dt;
    rho2=0.5*(rho+rho1);
    rho=rho-1i*(H*rho2-rho2*H)*dt+...
        (0.5*Gamma*(sigma*rho2*sigma'-sigma'*sigma*rho2+sigma*rho2*sigma'-rho2*sigma'*sigma))*dt;
end
%-------------------Ploting the result
plot(Time,real(n_s),'LineWidth',1.5)
xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
ylabel('$<\sigma^\dagger \sigma>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
ax = gca;
ax.FontSize = 20;
pbaspect([1 1 1])