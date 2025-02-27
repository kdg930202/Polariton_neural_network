clearvars
kk=1;
load("rho_s.mat")
%%%

which_rho = rho_s(:,:,5);
d=8;
dt = 0.05;
T = 0:dt:10-dt;
T = round(T,1);
time = 6;
tau = 1;
index_lmit=find((T>=time)&(T<time+tau));
J = 1;
trun_de = 0;
W = 1*rand(1,2);
pump = 3;
gamma = 1;
P = 0.1;
TD = 1;



d = 8; %dimension of the annihilation and creation operator
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


rho = kron(sparse(which_rho), rho_R);


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
                                         + 2*b2'*b2*rho3*b2'*b2 - rho3*b2'*b2*b2'*b2 - b2'*b2*b2'*b2*rho3);


    rho = rho + 1/6*dt*(K1+2*K2+2*K3+K4);

     if t>=index_lmit(1)
        t_after(kk)=T(t);
         rho_r_after_trace(:,:,kk)=partial_trace_source_out(d-1,2,1,2,full(rho))/trace(partial_trace_source_out(d-1,2,1,2,full(rho)));
          kk=kk+1;
     end 
end

%testing the rho_r_after
ss=size(rho_r_after_trace);
b1_after=kron(sig_m,I_b);index_lmit=find((T>=time)&(T<time+tau));
b2_after=kron(I_b,sig_m);index_lmit=find((T>=time)&(T<time+tau));
for ii=1:ss(3)
    n1_after(ii)=trace(b1_after'*b1_after*rho_r_after_trace(:,:,ii));
    n2_after(ii)=trace(b2_after'*b2_after*rho_r_after_trace(:,:,ii));
end

save('output.mat','b1_after','b2_after',' ')

figure()
plot(T, abs(n1))
hold on
plot(t_after,real(n1_after),'*')
%plot(T,abs(n2))