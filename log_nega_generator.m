% theta : 0~2pi
% s : 0.8~0.95
% phi : 0.5+-pi/10


% [v_min, N] = log_negativity(theta, s, phi)
clearvars
close all

num = 100;

theta_a = 0;
theta_b = 2*pi;
theta = (theta_b-theta_a).*rand(num,1) + theta_a;

s_a = 0.1;
s_b = 0.95;
s = (s_b-s_a).*rand(num,1) + s_a;


phi_l = ones(1,num)*(0.5+pi/10);
phi_r = ones(1,num)*(0.5-pi/10);
phi = [phi_l, phi_r];

for i = 1:num
    [V_min(i), NN(i)] = log_negativity(theta(i), s(i), phi(i));

end