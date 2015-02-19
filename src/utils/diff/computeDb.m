function [ Db ] = computeDb( myNEA , d , x )

[m,n] = size(x);
if n ~= 1
   error('In computing D(x), x should be a column vector')
end

q  = x(1:m/2,1);
dq = x(m/2+1:m,1);
myNEA  = myNEA.setState(q, dq);
D = myNEA.D.matrix;
b = myNEA.b.matrix;
Db = D * d + b;
end

