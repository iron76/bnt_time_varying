function [ Db ] = computeDb( myNEA , d , x )

[m,n] = size(x);
if n ~= 1
   error('In computing D(x), x should be a column vector')
end
NB = m/2;

q  = x(1:m/2,1);
dq = x(m/2+1:m,1);
myNEA  = myNEA.setState(q, dq);
if ~isa(myNEA, 'ANEA') && ~isa(myNEA, 'DANEA') 
   D = sparse(myNEA.iDs, myNEA.jDs, myNEA.Ds, 19*NB, 26*NB);
   b = sparse(myNEA.ibs, ones(size(myNEA.ibs)), myNEA.bs, 19*NB, 1);
   D = D(:, myNEA.id);
else
   D = myNEA.D.matrix;
   b = myNEA.b.matrix;
end

Db = D * d + b;
end

