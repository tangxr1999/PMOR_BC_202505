function F = first_lowrank_Q(A,N_j1,N_j2,B,C,alpha,N,g)
[~,u1] = size(B);
[~,u2] = size(C);
%u=2;
t = -2*alpha;
l_i = Laguerre(t,N);  %计算Laguerre函数l_i(t)也就是论文中的l_i(-2alpha)
f_i = cell(1,N);
L_i = cell(1,N);
X_i = cell(1,N);
H = cell(1,u1);
for k = 1:N
    f_i{1,k} = sqrt(2*alpha)*exp(alpha)*l_i{1,k};
end
%算出了fi

n1=size(A,1);
I=eye(n1);
L_i{1,1} = -A+1/sqrt(2*alpha)*I*f_i{1,1};
L_i{1,2} = A+1/sqrt(2*alpha)*I*(f_i{1,2}-2*f_i{1,1});
for h = 3:N
    L_i{1,h} = 1/sqrt(2*alpha)*I*(f_i{1,h}-2*f_i{1,h-1}+f_i{1,h-2});
end  %算出了Li

X_i{1,1} = sqrt(2*alpha)*L_i{1,1}^(-1);
X_i{1,2} = -L_i{1,1}^(-1)*L_i{1,2}*X_i{1,1};

for p = 3:N 
    L_ii= cell(1,p-1);
    X_ii = cell(p-1,1);
    for j = 1:p-1
    L_ii{1,j} = L_i{1,p+1-j};
    X_ii{j,1} = X_i{1,j};
    end
    L_iii = cell2mat(L_ii);
    X_iii = cell2mat(X_ii);
    XD = L_iii*X_iii;
    X_i{1,p} = -1*L_i{1,1}^(-1)*XD;    %用迭代计算X_i
end

% DD = M_integral(alpha,N);  %计算kronecker积那一块的
F = cell(1,g);  %F=[F1,F2];F1=FF1*D;FF1=[FF11,FF12,,,FF1N];F11=X_1{1,1}*B;
FF1 = cell(1,N);
FF2 = cell(1,N);

for d=1:g
    if d==1
        DD_1 = M1_integral(alpha,N,u2);%求Mi
        FF1{1,1} = X_i{1,1}*C;
        for m=2:N
            LL_i=cell(1,m-1);%N-1维行向量
            FF_i=cell(m-1,1);%N-1维列向量
            for mm=1:m-1
                LL_i{1,mm}=L_i{1,m+1-mm};
                FF_i{mm,1}=FF1{1,mm};
            end
            LL_ii = cell2mat(LL_i);
            FF_ii = cell2mat(FF_i);
            FD = LL_ii*FF_ii;
            FF1{1,m} = -L_i{1,1}^(-1)*FD;
        end
        FF1;
        F{1,d} = cell2mat(FF1)*DD_1;
    else
       
       % H = N_j1*F{1,d-1}
        H{1,1} = N_j1*F{1,d-1};
        H{1,2} = N_j2*F{1,d-1};
        FF2{1,1} = X_i{1,1}*cell2mat(H);%930*8
        DD_2 = M2_integral(alpha,N,d,u1);
        for z=2:N
            LLL_i=cell(1,z-1);
            FFF_i=cell(z-1,1);
            for zz=1:z-1
                LLL_i{1,zz}=L_i{1,z+1-zz};
                FFF_i{zz,1}=FF2{1,zz};
            end
            LLL_ii = cell2mat(LLL_i);%930*930
            FFF_ii = cell2mat(FFF_i);%930*8
            FDD=LLL_ii*FFF_ii;%930*8
            FF2{1,z}=-L_i{1,1}^(-1)*FDD;%930*8
        end
        FF2;
        F{1,d} = cell2mat(FF2)*DD_2;
    end
end
end