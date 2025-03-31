clearvars
close all

theta = 0:0.05:pi/2;
phi = 0:0.05:pi/2;


g2 = NaN(length(phi), length(theta));
ns = NaN(length(phi), length(theta));
n1_max = NaN(length(phi), length(theta));
% result_nu = zeros(length(phi), length(theta));
% n1_seleted = cell(length(phi), length(theta));
% n2_seleted = cell(length(phi), length(theta));

W = [0.5,0.5];
for i = 1:length(phi)
    disp(i)
    for j = 1:length(theta)

        [ns(i,j), g2(i,j), n1_selected{i,j}, n2_selected{i,j}] = g2_EQ13(W, phi(i), theta(j));
        % [ns(i,j), g2(i,j), n1_max(i,j)] = g2_only_n1_max(phi(i), theta(j));
        % g2(i,j) = g2_comparison_new_method(phi(i), theta(j));
        % [ns(i,j), g2(i,j)] = g2_ns_gs(phi(i), theta(j));

    end
end

%%
figure()
% subplot(1,3,1)
% pcolor(nth, alpha,result_an)
pcolor(theta,phi,abs(g2))
xlabel("theta")
ylabel("phi")
colorbar()
% title("Analytics")


% figure()
% subplot(1,3,1)
% pcolor(nth, alpha,result_an)
% pcolor(theta,phi,abs(n1_max))
% xlabel("theta")
% ylabel("phi")
% colorbar()

% figure()
% % subplot(1,3,1)
% % pcolor(nth, alpha,result_an)
% pcolor(theta,phi,abs(ns))
% xlabel("theta")
% ylabel("phi")
% colorbar()
% % title("Analytics")