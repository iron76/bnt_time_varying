function [a,m]=lu_cost(A)
% LU factorization computational cost.
% Input:
% A    the sparse matrix to be LU factorised
% Output:
% a    the number of additions
% m    the number of multiplications
%
% Francesco Nori, 06/21/16

a = 0;
m = 0;

[mA,nA] = size(A);
if mA~=nA
   error('The matrix to be LU-factorized must be square');
end

R = A;
P = cell(nA-1,1);
for k=1:nA-1
   if (R(k,k) == 0) error('Pivoting is needed!'); end
   r = k+1:nA;
   R(r,k) = R(r,k)/R(k,k);
   m = m + nnz(R(r,k));
   R(r,r) = R(r,r) - R(r,k)*R(k,r);
   a = a + nnz(R(r,k).*R(k,r)');
   m = m + nnz(R(r,k).*R(k,r)');
   P{k}   = R(r,r);
end
