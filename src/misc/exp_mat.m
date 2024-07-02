function y=exp_mat(A,x)
% EXP_MAT This function returns the exponential of matrix times a scalar x
% (step size in general) using precise integration scheme. More details
% about this integration scheme can be found at: W. Zhong and F.W. Williams
% A Precise Time Step Integration Method
%
% y = exp_mat(A,x)
%
% A: matrix
% x: step size
% y: output

N  = 20;
dx = x/2^N;
n  = size(A,1);

T1 = A*dx+(A*dx)^2*(eye(n)+A*dx/3+(A*dx)^2/12)/2;
for k=1:N
    T2=2*T1+T1*T1;
    T1=T2;
end
y = eye(n)+T1;

end