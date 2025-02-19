%Hamiltonian:
%H=(Delta1)*(a1^t*a1)+(Delta2)*(a2^t*a2)+
%hbar*J1*(a1^t*a2+a2^t*a1)+
%0.5*U1*(a1^t*a1^t*a1*a1)+
%hbar*(E1)*(a1^t+a1)+hbar*(E2)*(a2^t+a2)
clc;
clear;
%----constants
gamma1=1;% all the parametrs are scaled in gamma1
gamma2=0.1;
dDel = 0.01;
Delta1_vector=-5:dDel:5-dDel;
Delta2=1;
% E1_vector=0.01:0.02:0.5;% keep it in the low pumping regime
% E1 = omega_a, E2 = omega_b
E1_vector = 0.25;
phi=pi/2;
E2=0.4*exp(1i*phi);
J=1.6;%in the strong coupling regime it should be larger than decay rate
U=1;%
%------
Na1 = 7;
Na2 = 7;
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
%saving
save('rho_s.mat', 'rho_s')% It would be a huge file. It is 4D array. 
save('g2_a.mat','g2_a','Delta1_vector','E1_vector')
%After running code, you can plot g2 for different value of Deta2 or E1;

% figure;
% surf(Delta1_vector,E1_vector,(real(g2_a)));
% shading interp;
% view(0,90)
% colorbar
% xlabel('\Delta_a/\gamma_1');
% ylabel('\Omega_a/\gamma_1');
% title('g^{(2)}_a')
%str1=append('g^2_a(0) for U2=',U2_char);
% figure;
% surf(Delta1_vector,E1_vector,(real(g2_b)));
% shading interp;
% view(0,90)
% colorbar
% xlabel('E_1/\gamma_2')
% ylabel('\Delta_2/gamma_2');
% title('g^{(2)}_b')