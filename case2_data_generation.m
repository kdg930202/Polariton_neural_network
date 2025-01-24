clearvars
clc

r = 0.1:0.01:1;

for i=1:length(r)
    [g2_m1_psi_1(i)] = case2(r(i),1);
    % [g2_m3_psi_1(i)] = case2(r(i),3);
    % [g2_m5_psi_1(i)] = case2(r(i),5);
end


figure()
plot(r,g2_m1_psi_1)
hold on
% plot(r,g2_m3_psi_1)
% plot(r,g2_m5_psi_1)
% subplot(1,2,1)
% plot(alpha,g2_m1_psi_1)
% hold on
% plot(alpha,g2_m3_psi_1)
% plot(alpha,g2_m5_psi_1)
% legend(["m=1","m=3","m=5"])
% subplot(1,2,2)
% plot(alpha,g2_m1_psi_2)
% hold on
% plot(alpha,g2_m3_psi_2)
% plot(alpha,g2_m5_psi_2)
% legend(["m=1","m=3","m=5"])

