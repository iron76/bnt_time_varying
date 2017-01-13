function [ a, m ] = fw_cost( L, b )
% Forward substitution computational cost.
% Computes the computational cost of the forward 
% substitution to solve the linear system written as
%
%          L x = b
%
% Input:
% L    lower triangular matrix
% b    right-hand side of the linear system
% Output:
% a    the number of additions
% m    the number of multiplications
%
% Francesco Nori, 06/21/16

if norm(full(triu(L,1))) ~= 0
   error(['The provided L matrix in the fb_cost is not lower-triangular. Error is: ', num2str(norm(full(triu(L,1))))])
end


% Given a matrix L the solution of L x = b is recursively computed as:
%
%          xi = (bi - sum_{j=1}^{i-1} Lij xj)/Lii 
%
% with i = 1, 2, ..., N


a = nnz(b);
m = nnz(tril(L,-1));
a = a + nnz(tril(L,-1));
m = m + nnz(diag(L));
end

