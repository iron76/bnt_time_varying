function res = testSubmatrixSparse

res = 0;

A            = zeros(5, 10);
A1           = rand(2,3);
A2           = rand(1,7);
A3           = rand(2,7);
A(1:2,  1:3) = A1;
A(5  , 4:10) = A2;
A(3:4, 4:10) = A3;

As  = submatrix(      [2 2 1]', [3 7]', A);
Ass = submatrixSparse([2 2 1]', [3 7]', [1 3 2]', [1 2 2]');

Ass = set(Ass, A1, 1, 1);
Ass = set(Ass, A2, 3, 2);
Ass = set(Ass, A3, 2, 2);

for i = 1 : 3
   for j = 1 : 2
      if(norm(Ass(i,j)-As(i,j)) ~=0 )
         disp('[SUBMATRIX] Something wrong with submatrix access')
         res = 1;
      end
   end
end

if(norm(Ass.matrix-As.matrix) ~=0 )
   disp('[SUBMATRIX] Something wrong with the whole access')
   res = 1;
end


% As = set(As, A, [1 2 3], [1 2]);
% 
% res = 0;
% 
% [I, J] = indeces(As, 1, 1);
% if (norm(A(I,J)-A(1:2, 1:3))) || (norm(As(1,1)-A(I,J)))
%    disp('Something wrong with getting A(1,1)');
%    res = 1;
% end
% 
% [I, J] = indeces(As, 1, 2);
% if (norm(A(I,J)-A(1:2, 4:10))) || (norm(As(1,2)-A(I,J)))
%    disp('Something wrong with getting A(1,2)');
%    res = 1;
% end
% 
% [I, J] = indeces(As, 2, 1);
% if (norm(A(I,J)-A(3:4, 1:3))) || (norm(As(2,1)-A(I,J)))
%    disp('Something wrong with getting A(1,1)');
%    res = 1;
% end
% 
% [I, J] = indeces(As, 2, 2);
% if (norm(A(I,J)-A(3:4, 4:10))) || (norm(As(2,2)-A(I,J)))
%    disp('Something wrong with getting A(1,2)');
%    res = 1;
% end
% 
% [I, J] = indeces(As, 3, 1);
% if (norm(A(I,J)-A(5:5, 1:3))) || (norm(As(3,1)-A(I,J)))
%    disp('Something wrong with getting A(1,1)');
%    res = 1;
% end
% 
% [I, J] = indeces(As, 3, 2);
% if (norm(A(I,J)-A(5:5, 4:10))) || (norm(As(3,2)-A(I,J)))
%    disp('Something wrong with getting A(1,2)');
%    res = 1;
% end
% 
% As = set(As, zeros(1,7),3,2);
% if norm(As(3,2)-zeros(1,7))
%    disp('Something wrong with setting A(3,2)');
%    res = 1;
% end
% 
% if norm(As(1:2,1)-[As(1,1); As(2,1)])
%    disp('Something wrong with reading A(1:2,1)');
%    res = 1;
% end
% 
% if norm(As(2,[1 2])-[As(2,1) As(2,2)])
%    disp('Something wrong with reading A(2,[1 2])');
%    res = 1;
% end
% 
% if norm(As([1 3],2)-[As(1,2); As(3,2)])
%    disp('Something wrong with reading A([1 3],2)');
%    res = 1;
% end
% 
% [a,b] = size(As([1 3], 2));
% B     = rand(a,b);
% As    = set(As, B, [1 3], 2);
% 
% if norm(As([1 3],2)-B)
%    disp('Something wrong with setting A([1 3],2)');
%    res = 1;
% end
% 
% Bs = As([3 1], [1 2]);
% if norm(As([3 1],[1 2])-Bs)
%    disp('Something wrong with copying A([3 1],[1 2])');
%    res = 1;
% end


