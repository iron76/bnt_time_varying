% submatrix access matrix sub-blocks with arbitrary dimension
%
% This class allows to access a given matrix A subdivided into blocks of
% arbitrary (and not necessarily constant dimension). Lets assume that we
% want to access a matrix A with 'm' rows and 'n' columns. Rows are
% subdivided into M blocks of dimensions m_1, ..., m_M. Columns are
% subdivide into blocks of dimensions n_1, ..., n_N. Of course, we have
% that m = m_1 + ... + m_M and n = n_1 + ... + n_N. The class can be
% instantiated with the following definition:
%
% As = submatrix(A, [m_1 ... m_M], [n_1 ... n_N])
%
% The submatrix A_ij corresponding to the blocks m_i, n_j can be obtained
% by simply using:
%
% A_ij = As(i,j)
%
% The indeces to access the matrix A_ij are obtained as follows:
%
% [I, J] = indeces(As, i, j);
%
% so that we have A(I,J) = A_ij.
%
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef submatrix
   properties
      A, m, n, cm, cn
   end
   methods
      function b = submatrix(A, m, n)
         b.A  = A;
         b.m  = m;
         b.n  = n;
         b.cm = [0, cumsum(b.m)];
         b.cn = [0, cumsum(b.n)];
         
         
         [mA, nA] = size(A);
         if mA ~= b.cm(end)
            error('when calling submatrix(A, m, n) the sum(m) should equal the number of rows in A')
         elseif nA ~= b.cn(end)
            error('when calling submatrix(A, m, n) the sum(n) should equal the number of columns in A')
         end
      end
      
      function B = subsref(obj,S)
         if strcmp(S.type, '()') && length(S.subs) == 2
            [I, J] = obj.indeces(S.subs{1}, S.subs{2});
            B = obj.A(I, J);
            %          elseif strcmp(S.type, '()') && length(S.subs) == 4
            %             B = obj.A(S.subs{3}, S.subs{4});
         else
            error('The sumatrix class can be accessed only as A(i,j).')
         end
      end
      
      function [I, J] = indeces(obj, i, j)
         if i < 0 || i > length(obj.m) || j < 0 || j > length(obj.n)
            error('Trying to access the sumatrix outside its definition')
         else
            I = 1 + obj.cm(i) : obj.cm(i+1);
            J = 1 + obj.cn(j) : obj.cn(j+1);
         end
      end
      
      function disp(b)
         for i = 1 : length(b.m)
            fprintf('\n A(%d,*): \n \n', i)
            [I,~] = b.indeces(i, 1);
            for k = 1 : length(I)
               for j = 1 : length(b.n)
                  [I,J] = b.indeces(i, j);
                  %for k = 1 : length(I)
                  for h = 1 : length(J)
                     fprintf('%s ', num2str(b.A(I(k),J(h)), '%1.4f  '))
                  end
                  fprintf(' | ')
               end
               fprintf('\n')
            end
         end
      end % disp
      
      
      function obj = set(obj, Aij, i, j)
         [I,J] = obj.indeces(i, j);
         [mAij,nAij] = size(Aij);
         if length(I)~=mAij || length(J)~=nAij
            error('when setting Aij on As its dimension should match the dimensions of I and J in [I, J] = indeces(As, i, j)')
         else
            obj.A(I,J) = Aij;
         end
      end
      
   end
end