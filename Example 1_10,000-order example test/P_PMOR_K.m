function [V,odr,Er,Ar,Br,Cr]=P_PMOR_K(n,m,h,r,E,A,B,C,U,pK,epsilon)
% Generates the reduced discrete-time parametric system [Er,Ar,Br,Cr] by
% the P_PMOR_C method propsoed in Algorithm 2. 
% INPUT: n:the order of the original discrete-time parametric system;
%        m:the order of input
%        h:the index of truncating the Taylor series of system matrix A(p)
%        r:the index of truncating the expansion of the variable
%          wj(k)(j=0,1,...,h-1) and input u(t) in the space spanned by the
%          Krawtchouk polynomials
%        [E,A,B,C]: the system matrices of the original discrete-time
%                   parametric system, where A is a cell stored the first h 
%                   Taylor coefficients of the system matrix A(p)
%        U: the first r expansion coefficients of the input u(t) in the 
%           space spanned by the Krawtchouk polynomials
%        pK: the Krawtchouk polynomials parameter
%        epsilon：the positive constant for Algorithm 2, 0<epsilon<<1
% OUTPUT:V: the projection matrix obtained by Algorithm 2
%        odr: the order of the reduced discrete-time parametric system
%        [Er,Ar,Br,Cr]: the system matrices of the reduced discrete-time parametric system
   barE=kron(eye(h),E);
    barB=zeros(n*h,m);
    barB(1:n,:)=B;
    G=kron(eye(h),A{1});
    for i=1:1:h-1
        G=G+kron(diag(ones(h-i,1),-i),A{i+1});
    end
    del=epsilon^(1/r);
    Beta=(pK).^([0 0:1:r-2]);
    D=diag(del.^(0:1:r-1));
    F=conj(fft(eye(r)));
    Bu=flip(barB*U(1:r,:)',2);
    V=zeros(n*h,r);
    parfor i=0:r-1
        gi=(F(i+1,:)*D*Beta')*barE-G;
        bi=(1/sqrt(r))*Bu*D*F(i+1,:)';
        V0(:,i+1)=gi\bi;
    end
    T=(1/sqrt(r))*inv(D)*conj(F);
    for i=1:1:r
        V(:,i)=sum(V0*diag(T(i,:)),2);
    end
    V=reshape(V,n,r*h);
    [V,~]=qr(V,0);
    odr=size(V,2);
    Er=V'*E*V;Br=V'*B;Cr=C*V;
    for i=1:h
        Ar{i}=V'*A{i}*V;
    end
end