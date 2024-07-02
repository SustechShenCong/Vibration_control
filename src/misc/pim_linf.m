function y = pim_linf(A,f,x0,dt,t0,tf)
% PIM_LINF Forward simulation of a LTV system with external forcing. We use
% linear interpolation to approximate the external forcing here. The
% exponential of matrix is obtained via the precise integration method.
%
% y = pim_linf(A,f,x0,dt,tf)
%
% We solve for \dot{x} = Ax+f(t), x(0)=x0, 0<=t<=t_f. Here f is a
% functional handle
%

% calculation of transition matrix
s = exp_mat(A,dt);

% forward simulation via iteration
n = (tf-t0)/dt+1; n = round(n);
m = length(x0);
y = zeros(m,n);
y(:,1) = x0;

for k=1:n-1
    f0 = f(t0+(k-1)*dt);
    f1 = (f(t0+k*dt)-f0)/dt;
    y(:,k+1) = s*(y(:,k)+A\(f0+A\f1))-A\(f0+A\f1+dt*f1);
end

end