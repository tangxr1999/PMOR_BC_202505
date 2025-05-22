function [V,odr,Er,Ar,Br,Cr]=P_PMOR_C(n,m,h,r,E,A,B,C,U,ac)
% Generates the reduced discrete-time parametric system [Er,Ar,Br,Cr] by
% the P_PMOR_C method propsoed in Algorithm 1. 
% INPUT: n:the order of the original discrete-time parametric system;
%        m:the order of input
%        h:the index of truncating the Taylor series of system matrix A(p)
%        r:the index of truncating the expansion of the variable
%          wj(k)(j=0,1,...,h-1) and input u(t) in the space spanned by the
%          Charlier polynomials
%        [E,A,B,C]: the system matrices of the original discrete-time
%                   parametric system, where A is a cell stored the first h 
%                   Taylor coefficients of the system matrix A(p)
%        U: the first r expansion coefficients of the input u(t) in the 
%           space spanned by the Charlier polynomials
%        ac: the Charlier polynomials parameter
% OUTPUT:V: the projection matrix obtained by Algorithm 1
%        odr: the order of the reduced discrete-time parametric system
%        [Er,Ar,Br,Cr]: the system matrices of the reduced discrete-time parametric system
    barE=kron(eye(h),E);
    barB=zeros(n*h,m);
    barB(1:n,:)=B;
    G=kron(eye(h),A{1});
    for i=1:1:h-1
        G=G+kron(diag(ones(h-i,1),-i),A{i+1});
    end
    Bu=barB*U(1:r,:)';
    V=zeros(n*h,r);
    EG=barE-G;
    parfor i=0:r-1  
        bi=Bu(:,r-i);
        for j=1:1:i
            ci=EG\Bu(:,r-i+j);
            di=barE*ci;
            for k=1:1:j-1
                ci=EG\ci;
                di=barE*ci;
            end
            bi=bi+(pochhammer(r-i,j)/(ac^j))*di;
        end
        V(:,i+1)=EG\bi;
    end
    V=reshape(V,n,r*h);
    [V,~]=qr(V,0);
    odr=size(V,2);
    Er=V'*E*V;Br=V'*B;Cr=C*V;
    for i=1:h
        Ar{i}=V'*A{i}*V;
    end
end

