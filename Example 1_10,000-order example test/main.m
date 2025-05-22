%%
clear all;
clc
load('DateODE1.mat')
% n = size(A0,1);
N = 2;
g = 2;
alpha = 0.0874;
tol = 1e-10;
m = size(B,2);
h = 4;
rr = 5;
ac = 0.05;
PK = 0.35; % the Krawtchouk polynomials parameter
epsilon = 1e-12; % the constant ϵ
%% Ap
p_value = 0.1;
e0 = 2;
Ap = @(p) A0+p*A1+p^2*A2+p^3*A3;
%% System discretization
dt = 1e-2;
Ed0 = @(p) speye(n)-dt*Ap(p)/2;
Ad0 = @(p) speye(n)+dt*Ap(p)/2;
Bd0 = dt*B;
Cd0 = C;

A = cell(1,h);
A{1,1} = speye(n,n)+dt*A0/2;
A{1,2} = dt*A1/2;
A{1,3} = dt*A2/2;
A{1,4} = dt*A3/2;
%%
Ad1 = @(p) (Ed0(e0)^(-1))*Ad0(p);
Bd1 = (Ed0(e0)^(-1))*Bd0;
Cd1 = Cd0;
%% Bilinearization of parametric systems
A_bl = (Ed0(e0)^(-1))*(speye(n)+dt*A0/2);
N1 = zeros(n);
N2 = dt*(Ed0(e0)^(-1))*A1/2;
N3 = dt*(Ed0(e0)^(-1))*A2/2;
N4 = dt*(Ed0(e0)^(-1))*A3/2;
B0 = zeros(n,3);
B_bl = [Bd1,B0];
C_bl = Cd1;
m_bl = size(B_bl,2);   %System Input
q_bl = size(C_bl,1);   %System Output
%%
tic
[Ap_r, Br, Cr,r] = pmorByBilinear(A_bl, N1, N2, N3, N4, B_bl, C_bl, Ad1, Bd1, Cd1, alpha, N, g, tol, p_value);
t_pmor=toc;
disp(['The reduced order system dimension is',num2str(r)])
%%
U = [uk_expansion(ac,rr)]';
tic
[V,odr,Er,Ar,Br_1,Cr_1] = P_PMOR_C(n,m,h,rr,Ed0(e0),A,Bd0,Cd0,U,ac);
t_PMOR_C=toc;
Aprr = @(p) Ar{1,1}+p*Ar{1,2}+p^2*Ar{1,3}+p^3*Ar{1,4};
Apr = Aprr(p_value);
Ap_2 = Er^(-1)*Apr;
Br_2 = Er^(-1)*Br_1;
%%
UU = [uk_expansion(epsilon,rr)]';
tic
[V11,odr11,Er11,Ar11,Br_11,Cr_11] = P_PMOR_K(n,m,h,rr,Ed0(e0),A,Bd0,Cd0,UU,PK,epsilon);
t_PMOR_K=toc;
Aprr11 = @(p) Ar11{1,1}+p*Ar11{1,2}+p^2*Ar11{1,3}+p^3*Ar11{1,4};
Apr1 = Aprr11(p_value);
Ap_22 = Er11^(-1)*Apr1;
Br_22 = Er11^(-1)*Br_11;
%%
y = solve(Ad1(p_value), Bd1, Cd1);
yr1 = solve(Ap_r, Br, Cr);
yr2 = solve(Ap_2, Br_2, Cr_1);
yr3 = solve(Ap_22, Br_22, Cr_11);
%%
T0=linspace(1,100,100);
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
T1=linspace(0,100,51);
T1(T1==0) = 1;
plot(T1,y(T1),'k* -',T1,yr1(T1),'ro -',T1,yr2(T1),'b -.',T1,yr3(T1),'gd--','markersize',7,'LineWidth',1.3)
legend('Original, order n=10000','PMOR-BL, order r=4','P-PMOR-C, order r=20','P-PMOR-K, order r=20')
title('Transient Response')
xlabel('Time(k)')
ylabel('y(k)')
%%
figure(2)
semilogy(T1,aerr1(T1),'ro -',T1,aerr2(T1),'b -.',T1,aerr3(T1),'gd--','markersize',7,'LineWidth',1.3)
ylim([1e-9,1e-4]);
legend('PMOR-BL, order r=4','P-PMOR-C, order r=20','P-PMOR-K, order r=20')
title('Absolute Errors')    
xlabel('Time(k)')
ylabel('|y(k)-yr(k)|')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T=linspace(1,1000,1000);
% for j=1:length(T)
%     aerr1(j) = norm(y(:,j) - yr(:,j));
%     aerr2(j) = norm(y(:,j) - yr_1(:,j));
%     err1(j) = aerr1(j)/abs(y(:,j));
%     err2(j) = aerr2(j)/abs(y(:,j));
% end
% yh=y(1,1:100);
% yrh=yr(1,1:100);
% yrh_1=yr_1(1,1:100);
% aerr1h=aerr1(1,1:100);
% aerr2h=aerr2(1,1:100);
% T=linspace(0,100,51);
% T(T==0) = 1;
% %%
% figure(1)
% plot(T,yh(T),'-o',T,yrh(T),'-s',T,yrh_1(T),'-p','markersize',8,'LineWidth',1.5)
% legend('Original, order n=10000','PMOR-BL, order r=4','P-PMOR-C, order r=20')
% title('Transient Response')
% xlabel('Time(k)')
% ylabel('y(k)')
% %%
% figure(2)
% semilogy(T,aerr1(T),'r-o',T,aerr2(T),'y-s')
% ylim([1e-9,1e-4]);
% legend('PMOR-BL, order r=4','P-PMOR-C, order r=20')
% title('Absolute Errors')    
% xlabel('Time(k)')
% ylabel('|y(k)-yr(k)|')