clearvars
close all

theta = 0:0.1:pi;
phi = 0:0.1:pi;


g2 = NaN(length(phi), length(theta));
% result_nu = zeros(length(phi), length(theta));
n1_seleted = cell(length(phi), length(theta));
n2_seleted = cell(length(phi), length(theta));


for i = 1:length(phi)
    for j = 1:length(theta)

        [g2(i,j), n1_selected{i,j}, n2_selected{i,j}] = g2_EQ13(phi(i), theta(j));
        % g2(i,j) = g2_comparison_new_method(phi(i), theta(j));

    end
end

%%
figure()
% subplot(1,3,1)
% pcolor(nth, alpha,result_an)
pcolor(theta,phi,g2)
xlabel("theta")
ylabel("phi")
colorbar()
% title("Analytics")