clearvars

num = 1000;
theta_a = 0;
theta_b = 2*pi;
theta = (theta_b-theta_a).*rand(num,1) + theta_a;

s_a = 0.7;
s_b = 0.9;
s = (s_b-s_a).*rand(num,1) + s_a;

% [n1,n2] = QRP_single_mode(theta, s)
for i=1:num
    display(i)
    [n1,n2] = QRP_single_mode(theta(i), s(i));
    ni(i) = n1;
    nf(i) = n2;
end
