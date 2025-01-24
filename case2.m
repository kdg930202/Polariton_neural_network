% clc
% clearvars
function [g2_pis_1] = case2(r,m)

d = 50; %dimension of the annihilation and creation operator
aa = diag(sqrt(1:d-1),1); %annihilation operator

I_a = eye(d);
a1 = aa;


% Matrix exponential : expm
S = expm(0.5*(r'*a1*a1 - r*a1'*a1')); 

adm = I_a;
% m = 0;
for i=1:m
    adm = a1*adm;
end
vacc = I_a(:,1);
psi_1 = adm*S*vacc;
psi_1 = psi_1/norm(psi_1);


% psi_2 = zeros(1,length(psi_1))';
% for n = 0:d
%     if 2*n+1+m >= d
%         break
%     end
% 
%     % I_a(:,2*n+1)
%     psi_2 = psi_2 + 1/sqrt(cosh(alpha))*(-exp(1i*theta)*tanh(alpha))^n * sqrt(factorial(2*n+m))/(2^n*factorial(n)) * I_a(:,2*n+1+m);
% end
% psi_2 = 1/sqrt(cosh(alpha))*psi_2;

% PSI = [psi_1,psi_2];

rho1 = psi_1*psi_1';
trace(rho1)
% rho2 = psi_2*psi_2';

g2_pis_1 = trace(a1'*a1'*a1*a1*rho1)/trace(a1'*a1*rho1)^2;
% g2_pis_2 = trace(a1'*a1'*a1*a1*rho2)/trace(a1'*a1*rho2)^2;

end
%%
% clc
% m=0;
% I_a(:,m+1)

% Bm_sum = 0;
% for n=1:d
%     if 2*(n-1)+m>=d
%         break
%     end
%     % 2*(n-1)+m  
%     Bm_sum = Bm_sum + (r/2)^4*factorial(2*(n-1)+m)/(factorial(n-1)^2);
%     [2*(n-1)+m,n-1,factorial(2*(n-1)+m)/(factorial(n-1)^2)]
% end
% Nm = 1/sqrt(Bm_sum);
% psi = zeros(1,d)';
% for n=1:d
%     if 2*(n-1)+m+1>=d
%         break
%     end
%     % 2*(n-1)+m
%     Bm = (r/2)^2*sqrt(factorial(2*(n-1)+m))/(factorial(n-1));
%     psi = psi + Nm*Bm*I_a(:,2*(n-1)+m+1);
%     % [2*(n-1)+m,2*(n-1)+m+1,Bm]
% end
%%