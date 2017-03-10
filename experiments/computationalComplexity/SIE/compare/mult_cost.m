function [ a, m ] = mult_cost( A, B)
% Matrix multiplication computational cost.
% Computes the computational cost of the mutiplication
% between two matrices:
%
%         A B = C
%
% Input:
% A    first matrix (left-hand side)
% B    second matrix (right-hand side)
% Output:
% a    the number of additions
% m    the number of multiplications
%
% Francesco Nori, 06/21/16

[mA,nA] = size(A);
[mB,nB] = size(B);

if nA ~= mB
   error('A should have as many columns as the number of B columns')
end

% Given the matrices A and B, the (i,j) element in C is:
%
%          Cij = A(i,:) * B(:,j)
%
% with i,j = 1, 2, ..., N

a = 0;
m = 0;
for i = 1 : mA
   for j = 1 : nB
      z = nnz(A(i,:)' .* B(:,j));
      if z>1
         a = a + z-1;
      elseif z > 0
         m = m + z;
      end
   end
end

