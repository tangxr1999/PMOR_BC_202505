function V=uk1_expansion(a,N)
%expansion coefficient for u(k)=exp(-k)
V=exp(-a)*exp(a/exp(1));
for k=1:N-1
    V=[V exp(-a)*(1-1/exp(1))*V(end)];
end
end