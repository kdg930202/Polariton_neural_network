clearvars

d = 10; %dimension of the annihilation and creation operator
aa = diag(sqrt(1:d-1),1); %annihilation operator
I_a = eye(d);
a1 = aa;

theta = 1.9605;
s = 0.1651;
phi = 0.5-pi/10;

alpha = s*sin(phi);
a = abs(alpha)*exp(1i*theta);

% Matrix exponential : expm
S = expm(0.5*(a'*a1*a1 - a*a1'*a1')); 

%%
n_bar = s^2 * (cos(phi))^2 ;
% n_bar = 0.5;

for i=1:d
    rho_th(i) = (1/(1 + n_bar)) * (n_bar/(1 + n_bar))^(i-1);
    % Rho_th = Rho_th + rho_th * kron(I_a(:,i+1),I_a(:,i+1));
end

rho_th = diag(rho_th);
[trace(rho_th),trace(rho_th*a1'*a1)] %should be 1 and n_bar

rho_sq_th = S*rho_th*S';
% rho_sq_th = S*rho_th*conj(S);
[trace(rho_sq_th),trace(rho_sq_th*a1'*a1)] %should be 1 and n_bar
