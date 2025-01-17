clc
clearvars

d = 6; %dimension of the annihilation and creation operator
aa = diag(sqrt(1:d-1),1); %annihilation operator
s = 1;
theta = pi/4;

I_a = eye(d);
a1 = aa;
phi = 0.5+pi/10;
alpha = s*sin(phi);
r = abs(alpha)*exp(1i*theta);

% Matrix exponential : expm
S = expm(0.5*(r'*a1*a1/2 - r*a1'*a1'/2)); 

adm = I_a;
m = 2;
for i=1:m
    adm = a1'*adm;
end
vacc = I_a(:,1);
psi_1 = adm*S*vacc;

%%
clc
m=0;
ns = 0:2:d;
% I_a(:,m+1)

Bm_sum = 0;
for n=1:d
    if 2*(n-1)+m>d
        break
    end
    2*(n-1)+m  
    Bm_sum = Bm_sum + factorial(2*(n-1)+m)/(factorial(n-1)^2);
end
Nm = 1/sqrt(Bm_sum);
psi = zeros(1,d)';
for n=1:d
    if 2*(n-1)+m+1>d
        break
    end
    % 2*(n-1)+m
    Bm = (r/2)^2*factorial(2*(n-1)+m)/(factorial(n-1)^2);
    psi = psi + Nm*Bm*I_a(:,2*(n-1)+m+1);
end
%%