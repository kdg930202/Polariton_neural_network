%%Time evolution of a two level system under coherent pumping
clc;
clear;
%close all;
%-------------------Constants
Gamma_a=0.5;
Gamma_c=1;
E_a=1;
E_c=1;
Omega=1;
alpha=0;
E=0.2;
%-----------------Operators
N=3;
a1=diag(sqrt(1:N),1);
Ia=eye(N+1);
b1=[0 0;1 0];%sigma_m
Ib=eye(2);
a=kron(a1,Ib);
b=kron(Ia,b1);
bz=b'*b-b*b';
%------------------Hamiltonian
H=E_a*bz+E_c*a'*a+Omega*(a*b'+a'*b)+...
    alpha*a'*a'*a*a+E*(a'+a);
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
rho_c=p1*p1';

q1=zeros(2,1);q1(1)=1;
rho_a=q1*q1';

rho=kron(rho_c,rho_a);
%----------------------Dynamics
n_c=zeros(1,Nt);
n_a_0=zeros(1,Nt);
n_a_1=zeros(1,Nt);
for i=1:Nt
    %second order
    n_c(i)=trace(a'*a*rho);
    n_a_0(i)=trace(b*b'*rho);
    n_a_1(i)=trace(b'*b*rho);
    k1=-1i*(H*rho-rho*H)*h+...
        h*(0.5*Gamma_c*(a*rho*a'-a'*a*rho+a*rho*a'-rho*a'*a))...
        +h*(0.5*Gamma_a*(b*rho*b'-b'*b*rho+b*rho*b'-rho*b'*b));
    rho1=rho+0.5*k1;
    k2=-1i*(H*rho1-rho1*H)*h+...
        h*(0.5*Gamma_c*(a*rho1*a'-a'*a*rho1+a*rho1*a'-rho1*a'*a))...
        +h*(0.5*Gamma_a*(b*rho1*b'-b'*b*rho1+b*rho1*b'-rho1*b'*b));
    rho=rho+k2;
    disp(i)
end
%-------------------Ploting the result
figure;
plot(Time,real(n_c),'LineWidth',1.5);
hold on;
plot(Time,real(n_a_0),'LineWidth',1.5);
plot(Time,real(n_a_1),'LineWidth',1.5);
xlabel('Time','FontSize', 22,'FontName', 'Times New Roman')
ylabel('$<a^\dagger a>$','Interpreter', 'latex', 'FontSize', 22,'FontName', 'Times New Roman')
ax = gca;
ax.FontSize = 20;
pbaspect([1 1 1])