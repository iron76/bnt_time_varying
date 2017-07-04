function [L,U,P]=my_lu(A)
% LU factorization of an n by n matrix A
% using outher product LU
% A is factored as A = L*U
% Output:
% L is lower triangular with the main diagonal part = 1s.
% U is upper triangular and is stored in the original mtx A
% and must be zeroed out to get U
% Francesco Nori, 23 Jun 2016

[m,n] = size(A);
if m~=n
   error('The matrix to be LU-factorized must be square');
end
R = A;
P = cell(n-1,1);
for k=1:n-1
   if (R(k,k) == 0) error('Pivoting is needed!'); end
   r = k+1:n;
   R(r,k) = R(r,k)/R(k,k);
   R(r,r) = R(r,r) - R(r,k)*R(k,r);
   P{k}   = R(r,r);
end
U = triu(R);
L = tril(R,-1) + eye(size(R));