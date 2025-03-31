clearvars
close all

% r = linspace(0,1,1000);
r=0;
% result_nu = zeros(length(phi), length(theta));
n_all = cell(1,length(r));
% n2_selected = cell(1,length(r));


for i = 1:length(r)
    i
    % [g2_m(i), n1_selected{i}, n2_selected{i}] = case2_data_gen(r(i));
    [n_all{i}, g2_m(i)] = case2_data_gen(r(i));
end

%%
% figure()
% subplot(1,3,1)
% pcolor(nth, alpha,result_an)
% pcolor(theta,phi,g2)
% xlabel("theta")
% ylabel("phi")
% colorbar()
% title("Analytics")