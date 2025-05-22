function [Ap_r, Br, Cr, r] = pmorByBilinear(A_bl, N1, N2, B_bl, C_bl, Ad1, B, C, alpha, N, g, tol, p_value)
FF = first_lowrank(A_bl, N1, N2, B_bl, alpha, N, g);
F=cell2mat(FF);
GG = first_lowrank_Q(A_bl', N1', N2', B_bl,  C_bl', alpha, N, g);
G=cell2mat(GG);
%%
[U, S, V] = svd((G'*F), 'econ');
rN = sprank(S);
Sigma = diag(S(1:rN,1:rN));
r_adp = 1;
while(2*sum(Sigma(r_adp+1:rN)) > tol)
    r_adp = r_adp + 1;
end
r = r_adp;
% r = 10;
Sr = S(1:r, 1:r);
Ur = U(:, 1:r);
Vr = V(:, 1:r);
T=F*Vr*Sr^(-1/2);
S=Sr^(-1/2)*Ur'*G';
%% ½µ½×ÏµÍ³
Ap=Ad1(p_value);
Ap_r=S*Ap*T;
Br=S*B;
Cr=C*T;
end