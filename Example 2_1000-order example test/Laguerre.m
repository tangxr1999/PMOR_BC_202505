function l_i = Laguerre(t,N)
l_i = cell(1,N);
l_i{1,1} = 1;
l_i{1,2} = 1-t;
for k = 2:N-1
    l_i{1,k+1} = (1+2*(k-1)-t)*l_i{1,k}-(k-1)^2*l_i{1,k-1};   %Laguerre多项式的递推公式
end
end