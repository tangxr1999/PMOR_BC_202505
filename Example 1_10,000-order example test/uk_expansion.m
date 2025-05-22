function V=uk_expansion(a,N)
%expansion coefficient for u(k)=cos(k)
V1=exp(-a)*exp(a*cos(1))*cos(a*sin(1));
V2=exp(-a)*exp(a*cos(1))*sin(a*sin(1));
for k=1:N-1
    V1=[V1 exp(-a)*((1-cos(1))*V1(end)+sin(1)*V2(end))];
    V2=[V2 exp(-a)*((1-cos(1))*V2(end)-sin(1)*V1(end))];
end
V=V2;
%%for u=sin(k) set 
%V=V2;
end