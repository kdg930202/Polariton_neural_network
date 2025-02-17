clc;
clear;
close all;
%-------------------RF
Delta_s=1;
gamma_s=1;
Omega_s=0.4;
%------------------Laser
Delta_a_vector=-3:0.02:3;
gamma_a=1;
Omega_a=0.3;
%------------
theta=pi/4;

for i=1:length(Delta_a_vector)
    Delta_a=Delta_a_vector(i);
    [g2_O1,g2_O2,rho_out,alpha_num,alpha_ana,pop_s,pop_a,n_O1,n_O2]=g2_mix(Delta_s,gamma_s,Omega_s,Delta_a,gamma_a,Omega_a,theta);
    g2_O1_4saving(i)=g2_O1;
    g2_O2_4saving(i)=g2_O2;
    rho_out_4saving(:,:,i)=rho_out;
    
    %density in laser |alpha|^2
    pop_laser_4saving(i)=pop_a;
    %density in RF
    pop_RF_4saving(i)=pop_s;
    %density of O1 out put
    pop_O1_4saving(i)=n_O1;
    %density of O2 out put
    pop_O2_4saving(i)=n_O2;
    
    %to check if N is selected correctly
    test_N(i)=abs(alpha_num)^2-abs(alpha_ana)^2;
end

semilogy(Delta_a_vector,real(g2_O1_4saving),'linewidth',2);
hold on
semilogy(Delta_a_vector,real(g2_O2_4saving),'linewidth',2);

xlabel('\Delat_a/ \gamma_a')
ylabel('g^{(2)}_{O_2},g^{(2)}_{O_1} ')
legend('g^{(2)}_{O_1}','g^{(2)}_{O_2}')

 figure;
 plot(Delta_a_vector,real(pop_O1_4saving),'linewidth',2);
 hold on
plot(Delta_a_vector,real(pop_O2_4saving),'linewidth',2);
xlabel('\Delat_a/\gamma_a')
ylabel('n_{O_2},n_{O_1} ')
legend('n_{O_1}','n_{O_2}')

figure;
%we select N to make this almost zero
plot(Delta_a_vector,test_N)
% hold on;
% plot(Delta_a_vector,pop_O1_4saving+pop_O2_4saving)