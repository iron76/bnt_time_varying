function [a,m] = chol_cost(A)
% Computational cost of the Cholesky factorization of 
% an n by n matrix A using the Gaxpy-Rich algorithm
% A is factored as A = L*L'
% Output:
% a: number of additions
% m: number of mutiplications 
% Francesco Nori, 23 Jun 2016

[m,n] = size(A);
if m~=n
   error('The matrix to be Cholesky-factorized must be square');
end
a = 0;
m = 0;
for j=1:n
   if j > 1
      A(j:n,j) = A(j:n,j) - A(j:n,1:(j-1))*A(j,1:(j-1))';
      for k = j : n
         if nnz(A(k,j))~=0 || nnz(A(k,1:(j-1))*A(j,1:(j-1))') ~=0 
            a = a + nnz(A(k,j)) + nnz(A(k,1:(j-1))*A(j,1:(j-1))') - 1 ;
         end
         m = m + nnz(A(k,1:(j-1))*A(j,1:(j-1))');
      end
   end
   A(j:n,j) = A(j:n,j)./sqrt(A(j,j));
   m = m + 1;
end
