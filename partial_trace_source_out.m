function rho_r=partial_trace_source_out(N_s,n_mode_ins,N_r,n_mode_inr,rho_total)
%-------------
% For the bosonic source mode, this number shows maximum occupation number is 7
%N_s=3;
% This is the dimension of source Hilbert space; We have two modes in
% source, and then it needs to be squared;
%n_mode_ins=2;
d_s=(N_s+1)^n_mode_ins;%dimension of Hilbert space in Bosonic source
%N_s + 1 = d
%
%For the fermionic reservoir mode, the occupation number is 2
%N_r=2;
%For the fermionic reservoir, this number is the number of the modes
%n_mode_inr=2;
d_r=(N_r+1)^n_mode_inr;%Dimension of the Hilbert space in the Fermionic Reservoir.
I_r=eye(d_r);
%--------------------------ket states of the Bosonic space |n1,n2>, the
%maximum of n1 and n2 is N_s; |n1,n2> ---> kron(ket_s_only1mode,ket_s_only1mode)
ket_0=zeros(N_s+1,1);
ket_0(1)=1;
for i=1:N_s+1
    ket_s_only1mode(:,i)=ket_0;
    ket_0(i+1)=1;
    ket_0(i)=0;
    if i==N_s+1
        break
    end
end

kk=1;
for i=1:N_s+1
    for j=1:N_s+1
        ket_n1n2(:,kk)=kron(ket_s_only1mode(:,i),ket_s_only1mode(:,j));
        kk=kk+1;
    end
end
%---------------------------------The partial trace
%Give the density matrix rho, we should calculate this summation:
%rho_r=sum_{n1,n2}<n1,n2|rho|n1,n2>
rho_r=zeros(d_r,d_r);

for i=1:d_s
        ket_n1n2_s=ket_n1n2(:,i);
        rho_r=rho_r+kron(ket_n1n2_s,I_r)'*rho_total*kron(ket_n1n2_s,I_r);
end

% kk=1;
% for i=1:N_s+1
%         ket_n1n2=ket_s_only1mode(:,i);
%         rho_r=rho_r+kron(I_r,ket_n1n2)'*rho_total*kron(I_r,ket_n1n2);
%         kk=kk+1;
% end
%disp(1)