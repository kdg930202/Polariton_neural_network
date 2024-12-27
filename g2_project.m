clearvars

d = 10; %dimension of the annihilation and creation operator
a = diag(sqrt(1:d-1),1); %annihilation operator
I = eye(d);
alpha = 0.5;

%%
% Coherent state
rho_coh = 0; %Initialization

for i = 0:d-1
    for j = 0:d-1
    rho_coh = rho_coh + exp(-alpha^2)*(alpha^i)*(conj(alpha)^j)/sqrt(prod(1:i)*prod(1:j))*I(:,i+1)*I(:,j+1)';
    end
end

trace(rho_Coh)

%% 
% Thermal state
nth = 0.5;
rho_th = 0;

for i=0:d-1
    rho_th = rho_th+(nth^i/(1+nth)^(i+1))*I(:,i+1)*I(:,i+1)';
end
trace(rho_th)

%%
rho_fock = I(:,d)*I(:,d)';