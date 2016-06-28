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

U = A;
L = eye(nA);
for k=1:nA
    if (U(k,k) == 0) Error('Pivoting is needed!'); end
    L(k+1:nA,k)=U(k+1:nA,k)/U(k,k);
    m = m + nnz(U(k+1:nA,k));
    for j=k+1:nA
        U(j,:)=U(j,:)-L(j,k)*U(k,:);
        if L(j,k) ~= 0 
           a = a + nnz(U(k,:));
           m = m + nnz(U(k,:));
        end
    end
end
U = triu(U);