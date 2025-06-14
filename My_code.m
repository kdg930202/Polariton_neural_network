clearvars
clc
%%
dim = 2;;
Is = eye(dim);
v1 = Is(:,1);
v2 = Is(:,2);

sp = v1*v2';
sm = v2*v1';


ome = 1;
E = 1;
gamma = 0.1;

H = ome*(sp + sm) + E*(v2*v2');
L_H = -1i*kron(Is,H) + 1i*kron(H.',Is);

L_I = gamma*( kron(conj(sm),sm) - 0.5*kron(Is,sm'*sm) - 0.5*kron(sm.'*conj(sm),Is));

L_total = L_H + L_I;
TOL = 1e-6;
[R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L_total);

psi_0 = v1; %excited state
rho_0 = psi_0 * psi_0'; %Initial density matrix
Nt = 1000;
ti = 0;
tf = 10;
dt = (tf-ti)/(Nt-1);
t = ti:dt:tf;

for n=1:length(t)
    rho = zeros(dim,dim);
    for i = 1:length(lambda_sort)
        Lk = L_sort{i};
        Rk = R_sort{i};
        ck = trace(rho_0*Lk);
        rho = rho + ck*exp(lambda_sort(i)*t(n))*Rk;
    end
    rho00(n) = rho(1,1);
    rho11(n) = rho(2,2); 
end

figure()
plot(t, abs(rho00))
hold on
plot(t, abs(rho11))
