% theta : 0~2pi
% s : 0.8~0.95
% phi : 0.5+-pi/10


% clearvars

function [v_min, N] = log_negativity(theta, s, phi)

d = 6; %dimension of the annihilation and creation operator
aa = diag(sqrt(1:d-1),1); %annihilation operator
I_a = eye(d);
a1 = aa;
a2 = aa;

A1 = kron(a1, I_a);
A2 = kron(I_a, a2);


% n1 = d-1;
% n2 = d-1;

% 
% s = 0.8;
% phi = 0.5+pi/10;
% theta = pi;
alpha = s*sin(phi);
a = abs(alpha)*exp(1i*theta);

% Matrix exponential : expm
S = expm(a*A1'*A2' - a'*A1*A2); 


n_bar = s^2 * (cos(phi))^2 ;

n11 = zeros(1,d^2);
n22_ini = 0:d-1;
n22 = n22_ini;

for i = 1:d
    n11(d*i-d+1:d*i) = i-1;
    
end

for i=1:d-1
    n22 = [n22,n22_ini];
end

for i=1:d^2
    rho_th(i) = (1/(1 + n_bar))^2 * (n_bar/(1 + n_bar))^(n11(i)+n22(i));
    % Rho_th = Rho_th + rho_th * kron(I_a(:,i+1),I_a(:,i+1));
end

rho_th = diag(rho_th);
trace(rho_th);

rho_sq_th = S*rho_th*S';
% rho_sq_th = S*rho_th*conj(S);
trace(rho_sq_th);




q1 = (A1 + A1')/sqrt(2);
q2 = (A2 + A2')/sqrt(2);
p1 = (A1 - A1')/(1i*sqrt(2));
p2 = (A2 - A2')/(1i*sqrt(2));

V11 = 1/2*trace(rho_sq_th*(q1*q1+q1*q1)) - trace(rho_sq_th*q1)*trace(rho_sq_th*q1);
V12 = 1/2*trace(rho_sq_th*(q1*p1+p1*q1)) - trace(rho_sq_th*p1)*trace(rho_sq_th*q1);
V13 = 1/2*trace(rho_sq_th*(q1*q2+q2*q1)) - trace(rho_sq_th*q2)*trace(rho_sq_th*q1);
V14 = 1/2*trace(rho_sq_th*(q1*p2+p2*q1)) - trace(rho_sq_th*p2)*trace(rho_sq_th*q1);

V21 = 1/2*trace(rho_sq_th*(p1*q1+q1*p1)) - trace(rho_sq_th*q1)*trace(rho_sq_th*p1);
V22 = 1/2*trace(rho_sq_th*(p1*p1+p1*p1)) - trace(rho_sq_th*p1)*trace(rho_sq_th*p1);
V23 = 1/2*trace(rho_sq_th*(p1*q2+q2*p1)) - trace(rho_sq_th*q2)*trace(rho_sq_th*p1);
V24 = 1/2*trace(rho_sq_th*(p1*p2+p2*p1)) - trace(rho_sq_th*p2)*trace(rho_sq_th*p1);

V31 = 1/2*trace(rho_sq_th*(q2*q1+q1*q2)) - trace(rho_sq_th*q1)*trace(rho_sq_th*q2);
V32 = 1/2*trace(rho_sq_th*(q2*p1+p1*q2)) - trace(rho_sq_th*p1)*trace(rho_sq_th*q2);
V33 = 1/2*trace(rho_sq_th*(q2*q2+q2*q2)) - trace(rho_sq_th*q2)*trace(rho_sq_th*q2);
V34 = 1/2*trace(rho_sq_th*(q2*p2+p2*q2)) - trace(rho_sq_th*p2)*trace(rho_sq_th*q2);

V41 = 1/2*trace(rho_sq_th*(p2*q1+q1*p2)) - trace(rho_sq_th*q1)*trace(rho_sq_th*p2);
V42 = 1/2*trace(rho_sq_th*(p2*p1+p1*p2)) - trace(rho_sq_th*p1)*trace(rho_sq_th*p2);
V43 = 1/2*trace(rho_sq_th*(p2*q2+q2*p2)) - trace(rho_sq_th*q2)*trace(rho_sq_th*p2);
V44 = 1/2*trace(rho_sq_th*(p2*p2+p2*p2)) - trace(rho_sq_th*p2)*trace(rho_sq_th*p2);

A = [V11, V12; V21, V22];
C = [V13, V14; V23, V24];
B = [V33, V34; V43, V44];

V = [A, C;C', B];

sig = det(A) + det(B) - 2*det(C);
v_min = sqrt(sig - sqrt(sig^2 - 4*det(V)))/sqrt(2);

N = -log2(2*v_min); % if N is positive, the state rho is entangled.

end