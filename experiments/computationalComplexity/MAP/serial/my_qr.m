function [ Q,R ] = my_qr( A )
%MY_QR QR factorisation of a matrix A (Givens algorithm)

[m,n] = size(A);
Q = eye(m);
R = A;

% |I 0 0 | |C11 C12 C13|   | C11  C12  C13|
% |0 A 0 | |C21 C22 C23| = |AC21 AC22 AC23| 
% |0 0 I | |C31 C32 C33|   | C31  C32  C33|
% 
% |C11 C12 C13| |I 0 0 |   | C11 AC12  C13|
% |C21 C22 C23| |0 A 0 | = | C21 AC22  C23| 
% |C31 C32 C33| |0 0 I |   | C31 AC32  C33|

for j = 1:n
   for i = m:-1:(j+1)
      [c,s] = givensrotation( R(i-1,j),R(i,j) );
      R([i-1, i], j:n) = [c -s; s c]'*R([i-1, i], j:n);
      Q(:, [i-1, i]) = Q(:, [i-1, i])*[c -s; s c];
   end
end

R = triu(R);
end


function [c,s] = givensrotation(a,b)
if b == 0
   c = 1;
   s = 0;
else
   if abs(b) > abs(a)
      r = a / b;
      s = 1 / sqrt(1 + r^2);
      c = s*r;
   else
      r = b / a;
      c = 1 / sqrt(1 + r^2);
      s = c*r;
   end
end

end