% alpha:Laguerre parameter
%     N:Number of Laguerre function expansion terms
%     u:Number of Inputs
function DD_1=M1_integral(alpha,N,u)
a = alpha;
syms s;
for k=1 : N
    for j=1 : N
    f(k,j)=(-1).^(j-k)*2*a*((a+s*1i).^(j-k-1)*(a-s*1i).^(k-j-1));
     H(k,j) = int(f(k,j), s, [-pi,pi]);
    G(k,j) =H(k,j)/(2*pi);
     M(k,j) = double(G(k,j));
    end 
end
MM=real(M);%取实
DD=chol(MM); %Cholesky分解
DD_1 = kron(DD,eye(u));%
end