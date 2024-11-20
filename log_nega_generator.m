% theta : 0~2pi
% s : 0.8~0.95
% phi : 0.5+-pi/10


% [v_min, N] = log_negativity(theta, s, phi)
clearvars
close all

num = 500;

theta_a = 0;
theta_b = 2*pi;
theta = (theta_b-theta_a).*rand(num,1) + theta_a;

s_a = 0.7;
s_b = 0.9;
s = (s_b-s_a).*rand(num,1) + s_a;


phi_l = ones(1,num/2)*(0.5+pi/10);
phi_r = ones(1,num/2)*(0.5-pi/10);
phi = [phi_l, phi_r];
half = 0.5*ones(1,num);

for i = 1:num
    [V_min(i), NN(i)] = log_negativity(theta(i), s(i), phi(i));

end


figure()
polarscatter(theta(1:num/2), real(V_min(1:num/2)'),'b')
hold on
polarscatter(theta(num/2+1:end), real(V_min(num/2+1:end)'),'r')
polarplot(half','k')
% [v_sample, NN_sample] = log_negativity(0, 0.5, 0.5+pi/10);

% polarscatter(, NN_sample)
% 

