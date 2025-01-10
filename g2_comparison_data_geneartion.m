clearvars
close all

theta = 0:0.1:pi;
phi = 0:0.1:2*pi;

result_an = zeros(length(phi), length(theta));
result_nu = zeros(length(phi), length(theta));

for i = 1:length(phi)
    for j = 1:length(theta)

        [result_an(i,j), result_nu(i,j)] = g2_comparison(phi(i), theta(j));
        g2(i,j) = g2_comparison_new_method(phi(i), theta(j));

    end
end

%%

figure()
subplot(1,3,1)
% pcolor(nth, alpha,result_an)
pcolor(theta,phi,result_an)
xlabel("theta")
ylabel("phi")
title("Analytics")
colorbar()
subplot(1,3,2)
% pcolor(nth, alpha,result_nu)
pcolor(theta,phi,result_nu)
colorbar()
xlabel("theta")
ylabel("phi")
title("Numerics")
subplot(1,3,3)
% pcolor(nth, alpha,result_nu)
pcolor(theta,phi,abs(g2))
colorbar()
xlabel("theta")
ylabel("phi")
title("New")