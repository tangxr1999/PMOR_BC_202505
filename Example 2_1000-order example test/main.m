%2018-A parameterized model order reduction method for parametric systems
%based on Laguerre polynomials Example 1
%%
clear all ;
clc
N = 2;
g = 2;
alpha = 0.2010;
tol = 1e-10;
rr = 4;
ac = 0.01; % the Charlier polynomials parameter
PK = 0.35; % the Krawtchouk polynomials parameter
epsilon = 1e-12; % the constant ϵ
%% A0,A1,B,C
k = 500;
n = 2*k;
a = linspace(-1e3,10,k);
b = linspace(10,1e3,k);
A0 = sparse(n,n);
for k1 = 1:k
    A0(2*k1-1:2*k1, 2*k1-1:2*k1) = [0, b(k1); -b(k1), 0];
end
A1 = sparse(n,n);
for k2 = 1:k
    A1(2*k2-1, 2*k2-1) = a(k2);
    A1(2*k2, 2*k2) = a(k2);
end
B = sparse(n,1);
for k3 =1:k
    B(2*k3-1, 1) = 2;
end
m = size(B,2);
C = sparse(1,n);
for k4 = 1:k
    C(1, 2*k4-1) = 1;
end
%% Ap
p_value = 0.5;
e0 = 2;
Ap = @(p) A0+p*A1;
h = 2; % Number of truncated items of Ap

%% System discretization
dt = 1e-2;
Ed0 = @(p) speye(n,n)-dt*Ap(p)/2;
Ad0 = @(p) speye(n,n)+dt*Ap(p)/2;
Bd0 = dt*B;
Cd0 = C;

AA = cell(1,h); %a cell stored the first h 
                %Taylor coefficients of the system matrix A(p)
AA{1,1} = speye(n,n)+dt*A0/2;
AA{1,2} = dt*A1/2;
%%
% E0 = speye(n,n)-dt*A0/2;
% Ad1 = @(p) (E0^(-1))*Ad0(p);
Ad1 = @(p) (Ed0(e0)^(-1))*Ad0(p);
Bd1 = (Ed0(e0)^(-1))*Bd0;
Cd1 = Cd0;
%% Bilinearization of parametric systems
A_bl = (Ed0(e0)^(-1))*(speye(n,n)+dt*A0/2);
N1 = zeros(n);
N2 = dt*(Ed0(e0)^(-1))*A1/2;
B0 = zeros(n,1);
B_bl = [Bd1,B0];
C_bl = Cd1;
m_bl = size(B_bl,2);   %System Input
q_bl = size(C_bl,1);   %System Output
%%
tic
[Ap_r, Br, Cr,r] = pmorByBilinear(A_bl, N1, N2, B_bl, C_bl, Ad1, Bd1, Cd1, alpha, N, g, tol, p_value);
t_pmor=toc;
disp(['The reduced order system dimension is',num2str(r)])
%%
U = [uk2_expansion(ac,rr)]';
tic
[V,odr,Er,Ar,Br_1,Cr_1]=P_PMOR_C(n,m,h,rr,Ed0(p_value),AA,Bd0,Cd0,U,ac);
t_pmor_c=toc;
Aprr = @(p) Ar{1,1}+p*Ar{1,2};
Apr = Aprr(p_value);
Ap_2 = Er^(-1)*Apr;
Br_2 = Er^(-1)*Br_1;
%%
UU = [uk2_expansion(epsilon,rr)]';
tic
[V11,odr11,Er11,Ar11,Br_11,Cr_11]=P_PMOR_K(n,m,h,rr,Ed0(p_value),AA,Bd0,Cd0,UU,PK,epsilon);
t_pmor_k=toc;
Aprr11 = @(p) Ar11{1,1}+p*Ar11{1,2};
Apr1 = Aprr11(p_value);
Ap_22 = Er11^(-1)*Apr1;
Br_22 = Er11^(-1)*Br_11;
%%
y = solve(Ad1(p_value), Bd1, Cd1);
yr1 = solve(Ap_r, Br, Cr);
yr2 = solve(Ap_2, Br_2, Cr_1);
yr3 = solve(Ap_22, Br_22, Cr_11);
%%
T0=linspace(1,500,500);
for j=1:length(T0)
    aerr1(j) = norm(y(:,j) - yr1(:,j));
    aerr2(j) = norm(y(:,j) - yr2(:,j));
    aerr3(j) = norm(y(:,j) - yr3(:,j));
    err1(j) = aerr1(j)/abs(y(:,j));
    err2(j) = aerr2(j)/abs(y(:,j));
    err3(j) = aerr3(j)/abs(y(:,j));
end
%%
figure(1)
T1=linspace(0,500,51);
T1(T1==0) = 1;
plot(T1,y(T1),'k* -',T1,yr1(T1),'ro -',T1,yr2(T1),'b -.',T1,yr3(T1),'gd--','markersize',7,'LineWidth',1.3)
legend('Original, order n=1000','PMOR-BL, order r=5','P-PMOR-C, order r=8','P-PMOR-K, order r=8')
title('Transient Response')
xlabel('Time(k)')
ylabel('y(k)')
%%
figure(2)
semilogy(T1,aerr1(T1),'ro -',T1,aerr2(T1),'b -.',T1,aerr3(T1),'gd--','markersize',7,'LineWidth',1.3)
ylim([1e-6,1e1]);
legend('PMOR-BL, order r=5','P-PMOR-C, order r=8','P-PMOR-K, order r=8')
title('Absolute Errors')    
xlabel('Time(k)')
ylabel('|y(k)-yr(k)|')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 原作图代码
% T=linspace(1,1000,1000);
% for j=1:length(T)
%     aerr1(j) = norm(y(:,j) - yr(:,j));
%     aerr2(j) = norm(y(:,j) - yr_1(:,j));
%     err1(j) = aerr1(j)/abs(y(:,j));
%     err2(j) = aerr2(j)/abs(y(:,j));
% end
% yh=y(1,1:500);
% yrh=yr(1,1:500);
% yrh_1=yr_1(1,1:500);
% aerr1h=aerr1(1,1:500);
% aerr2h=aerr2(1,1:500);
% T=linspace(0,500,51);
% T(T==0) = 1;
% %%
% figure(1)
% plot(T,yh(T),'k* -',T,yrh(T),'ro -',T,yrh_1(T),'b --','markersize',8,'LineWidth',1.0)
% legend('Original, order n=1000','PMOR-BL, order r=5','P-PMOR-C, order r=8')
% title('Transient Response')
% xlabel('Time(k)')
% ylabel('y(k)')
% %%
% figure(2)
% semilogy(T,aerr1(T),'ro -',T,aerr2(T),'b --')
% ylim([1e-6,1e1]);
% legend('PMOR-BL, order r=5','P-PMOR-C, order r=8')
% title('Absolute Errors')    
% xlabel('Time(k)')
% ylabel('|y(k)-yr(k)|')