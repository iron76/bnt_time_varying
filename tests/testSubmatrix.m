clear all
close all
clc

m = 5;
n = 10;
A  = rand(5, 10);
As = submatrix(A, [2 2 1], [3 7]);

[I, J] = indeces(As, 1, 1);
if (norm(A(I,J)-A(1:2, 1:3))) || (norm(As(1,1)-A(I,J)))
   error('Something wrong with getting A(1,1)');
end

[I, J] = indeces(As, 1, 2);
if (norm(A(I,J)-A(1:2, 4:10))) || (norm(As(1,2)-A(I,J)))
   error('Something wrong with getting A(1,2)');
end

[I, J] = indeces(As, 2, 1);
if (norm(A(I,J)-A(3:4, 1:3))) || (norm(As(2,1)-A(I,J)))
   error('Something wrong with getting A(1,1)');
end

[I, J] = indeces(As, 2, 2);
if (norm(A(I,J)-A(3:4, 4:10))) || (norm(As(2,2)-A(I,J)))
   error('Something wrong with getting A(1,2)');
end

[I, J] = indeces(As, 3, 1);
if (norm(A(I,J)-A(5:5, 1:3))) || (norm(As(3,1)-A(I,J)))
   error('Something wrong with getting A(1,1)');
end

[I, J] = indeces(As, 3, 2);
if (norm(A(I,J)-A(5:5, 4:10))) || (norm(As(3,2)-A(I,J)))
   error('Something wrong with getting A(1,2)');
end

As = set(As, zeros(1,7),3,2);
if norm(As(3,2)-zeros(1,7))
   error('Something wrong with setting A(3,2)');
end

if norm(As(1:2,1)-[As(1,1); As(2,1)])
   error('Something wrong with reading A(1:2,1)');
end

if norm(As(2,[1 2])-[As(2,1) As(2,2)])
   error('Something wrong with reading A(2,[1 2])');
end

if norm(As([1 3],2)-[As(1,2); As(3,2)])
   error('Something wrong with reading A([1 3],2)');
end

[a,b] = size(As([1 3], 2));
B     = rand(a,b);
As    = set(As, B, [1 3], 2);

if norm(As([1 3],2)-B)
   error('Something wrong with setting A([1 3],2)');
end

Bs = As([3 1], [1 2]);
if norm(As([3 1],[1 2])-Bs)
   error('Something wrong with copying A([3 1],[1 2])');
end


