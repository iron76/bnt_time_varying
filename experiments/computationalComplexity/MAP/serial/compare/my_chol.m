function L = my_chol(A)
% Cholesky factorization of an n by n matrix A
% using the Gaxpy-Rich algorithm
% A is factored as A = L*L'
% Output:
% L is lower triangular 
% Francesco Nori, 23 Jun 2016

[m,n] = size(A);
if m~=n
   error('The matrix to be Cholesky-factorized must be square');
end
for j=1:n
   if j > 1
      A(j:n,j) = A(j:n,j) - A(j:n,1:(j-1))*A(j,1:(j-1))';
   end
   A(j:n,j) = A(j:n,j)./sqrt(A(j,j));
end
L = tril(A);