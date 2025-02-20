function[O2,g2_O1,g2_O2,rho_out,alpha_num,alpha_ana,pop_s,pop_a,n_O1,n_O2]=g2_mix(Delta_s,gamma_s,Omega_s,Delta_a,gamma_a,Omega_a,theta)
%------------RF
b=[0 0;1 0];
Ib=[1 0;0 1];
H_s=(Delta_s)*(b'*b)+Omega_s*(b'+b);
L_s=(1/2)*gamma_s*(2*kron(conj(b),b)-kron(Ib,b'*b)-kron(b.'*conj(b),Ib));

LL=-1i*kron(Ib,H_s)+1i*kron(H_s,Ib)+L_s;

[rhoS_s,lambda0] = eigs(LL,1,'sm'); 
eigen0 = lambda0;
rhoS_s = reshape(rhoS_s,2,2);
rhoS_s = rhoS_s/trace(rhoS_s);
pop_s= trace(b'*b*rhoS_s);
g2_s=trace((b'*b*b'*b-b'*b)*rhoS_s)/(trace(b'*b*rhoS_s))^2;
%disp(g2_s);
%disp(pop_s);
%------------laser
Na =11;
dima = Na + 1;
a=spdiags(sqrt(0:Na)',1,dima,dima);
Ia = speye(dima);
%--------------------
H_a=(Delta_a)*(a'*a)+Omega_a*(a'+a);
L_a=(1/2)*gamma_a*(2*kron(conj(a),a)-kron(Ia,a'*a)-kron(a.'*conj(a),Ia));

L=-1i*kron(Ia,H_a)+1i*kron(H_a,Ia)+L_a;

[rhoS1,lambda0] = eigs(L,1,'sm'); 
eigen0 = lambda0;
rhoS1 = reshape(rhoS1,dima,dima);
rhoS1 = rhoS1/trace(rhoS1);
pop_a= trace(a'*a*rhoS1);
alpha_num=trace(a*rhoS1);
alpha_ana=Omega_a/(1i*(gamma_a/2)-Delta_a);
g2_a=trace((a'*a*a'*a-a'*a)*rhoS1)/(trace(a'*a*rhoS1))^2;
% disp(alpha_ana)
% disp(alpha_num)
%------------------------rho_out
rho_in=kron(rhoS1,rhoS_s);
%trace(rho_in)
%theta=pi/2;

U=expm(-1i*theta*(kron(a',b)+kron(a,b')));

%ss=U*U';

rho_out=U*rho_in*U';
%rho_out=rho_out/trace(rho_out);
% disp(trace(rho_out));
O2=1i*sin(theta)*kron(a,Ib)+cos(theta)*kron(Ia,b);
n_O2=trace(O2'*O2*rho_out);
g2_O2=trace((O2'*O2'*O2*O2)*rho_out)/(trace(O2'*O2*rho_out))^2;

O1=cos(theta)*kron(a,Ib)+1i*sin(theta)*kron(Ia,b);
n_O1=trace(O1'*O1*rho_out);
g2_O1=trace((O1'*O1'*O1*O1)*rho_out)/(trace(O1'*O1*rho_out))^2;
end
