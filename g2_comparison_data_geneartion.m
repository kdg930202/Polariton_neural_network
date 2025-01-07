clearvars
close all

nth = 0:0.01:3;
alpha = 0:0.01:3;

result_an = zeros(length(nth), length(alpha));
result_nu = zeros(length(nth), length(alpha));

for i = 1:length(nth)
    for j = 1:length(alpha)

        [result_an(i,j), result_nu(i,j)] = g2_comparison(nth(i), alpha(j));
    end
end

%%

figure()
subplot(1,2,1)
% pcolor(nth, alpha,result_an)
pcolor(result_an)
xlabel("nth")
ylabel("alpha")
title("Analytics")
colorbar()
subplot(1,2,2)
% pcolor(nth, alpha,result_nu)
pcolor(result_nu)
colorbar()
xlabel("nth")
ylabel("alpha")
title("Numerics")