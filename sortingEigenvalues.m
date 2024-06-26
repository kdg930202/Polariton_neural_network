function [R_sort,L_sort,lambda_sort] = sortingEigenvalues(dim,TOL,L) 

[R,DR] = eig(L);
[L,DL] = eig(L');
eig_R = diag(DR);
eig_L = diag(DL);
ind_RL = zeros(dim*dim ,2);
count = 1;

for n=1:dim*dim
    an = eig_R(n);
    for m=1:dim*dim
        bm = eig_L(m);
        if (abs(real(an)-real(bm))<TOL && abs(imag(an)-imag(bm))<TOL && count<=dim*dim)
            ind_RL(count ,1) = n;
            ind_RL(count ,2) = m;
            count = count + 1;
        end
    end
end
eig_L = eig_L(ind_RL(:,2)');
eig_R = eig_R(ind_RL(:,1)');
L = L(:,ind_RL(:,2)');
R = R(:,ind_RL(:,1)');
lambda = eig_R;
[~,ind] = sort(lambda);
lambda_sort = lambda(ind);
L = L(:,ind);
R = R(:,ind);

% Final R_sort and L_sort matrices
R_sort = cell(1,length(lambda_sort));
L_sort = cell(1,length(lambda_sort));
for k=1:length(lambda_sort)
    R_sort{k} = reshape(R(:,k),dim,dim);
    L_sort{k} = reshape(L(:,k),dim,dim);
    Rk = R_sort{k};
    Lk = L_sort{k};
    Ck = trace(Lk*Rk);
    Lk = Lk/sqrt(Ck);
    Rk = Rk/sqrt(Ck);
    R_sort{k} = Rk;
    L_sort{k} = Lk;
end
end