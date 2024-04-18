%Time evolution of a two level system under coherent pumping
clc;
clear;
%close all;
%-------------------Constants
Gamma=1;
Omega=1;
alpha=2;
E=0.2;
%-----------------Operators
N=4;
a=diag(sqrt(1:N),1);
Ia=eye(N+1);
%------------------Hamiltonian
H=Omega*a'*a+E*(a'+a)+alpha*a'*a'*a*a;
%-------------------Time
T_max=100;nkk=1;
h=0.01;Time=0:h:T_max;Nt=length(Time);
%-----------------Initial conditions
ket_gen=zeros(N+1,1);
p1=ket_gen;
p1(1)=1;
% p2=ket_gen;
% p2(2)=1;
% p3=ket_gen;
% p3(3)=1;
% p4=ket_gen;
% p4(4)=1;
 p5=ket_gen;
 p5(5)=1;
rho=p1*p1';
%----------------------Dynamics
n_s=zeros(1,Nt);
for i=1:Nt
    %second order
    n_s(i)=trace(a'*a*rho);
    k1=-1i*(H*rho-rho*H)*h+...
        h*(0.5*Gamma*(a*rho*a'-a'*a*rho+a*rho*a'-rho*a'*a));
    rho1=rho+0.5*k1;
    k2=-1i*(H*rho1-rho1*H)*h+...
        h*(0.5*Gamma*(a*rho1*a'-a'*a*rho1+a*rho1*a'-rho1*a'*a));
    rho=rho+k2;

end
%-------------------Ploting the result
figure;
plot(Time,real(n_s),'LineWidth',1.5)
xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
ylabel('$<a^\dagger a>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
ax = gca;
ax.FontSize = 20;
pbaspect([1 1 1])