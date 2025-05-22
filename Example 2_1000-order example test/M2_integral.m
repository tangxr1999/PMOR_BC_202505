% alpha:Laguerre parameter
%     N:Number of Laguerre function expansion terms
%     d:Item d of the Gram matrix
%     u:Number of Inputs
function DD_2=M2_integral(alpha,N,u,v)
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
DD_2 = kron(DD,eye(u*v));
end