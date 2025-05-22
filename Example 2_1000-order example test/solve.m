function y = solve(A, B, C)
n = size(A,1);
x0 = zeros(n,1);
maxstep = 500;
x_initial = x0;
y(1) = C*x_initial;
 for k = 2:maxstep
     % uu = [exp(-k*0.02)*sin(k*0.02)]';
     uu = cos(k);
     % uu = exp(-k);
     % uu = [exp(-k*0.02)*sin(k*0.02)]';
     x_initial(:,k) = A*x_initial(:,k-1)+B*uu;
     y(k) = C*x_initial(:,k);
 end
end