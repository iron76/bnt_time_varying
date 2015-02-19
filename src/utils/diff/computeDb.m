function [ Db ] = computeDb( mySNEA , d , x )

[m,n] = size(x);
if n ~= 1
   error('In computing D(x), x should be a column vector')
end

q  = x(1:m/2,1);
dq = x(m/2+1:m,1);
mySNEA  = mySNEA.setState(q, dq);
D = mySNEA.D.matrix;
b = mySNEA.b.matrix;
Db = D * d + b;
end

