clearvars

alpha = 0:0.01:1;

for i=1:length(alpha)
    [g2_m1_psi_1(i),g2_m1_psi_2(i)] = case2(alpha(i),1);
    [g2_m3_psi_1(i),g2_m3_psi_2(i)] = case2(alpha(i),3);
    [g2_m5_psi_1(i),g2_m5_psi_2(i)] = case2(alpha(i),5);
end


figure()
subplot(1,2,1)
plot(alpha,g2_m1_psi_1)
hold on
plot(alpha,g2_m3_psi_1)
plot(alpha,g2_m5_psi_1)
legend(["m=1","m=3","m=5"])
subplot(1,2,2)
plot(alpha,g2_m1_psi_2)
hold on
plot(alpha,g2_m3_psi_2)
plot(alpha,g2_m5_psi_2)
legend(["m=1","m=3","m=5"])

