function [ D ] = computeD( myNEA , q)

[~,n] = size(q);
if n ~= 1
   error('In computing D(x), x should be a column vector')
end

dq     = zeros(size(q));
myNEA  = myNEA.setState(q, dq);
if isa(myNEA, 'DANEA') || isa(myNEA, 'ANEA')
   D = myNEA.D.matrix;
else
   D = sparse(myNEA.iDs, myNEA.jDs, myNEA.Ds, 19*myNEA.IDmodel.modelParams.NB, 26*myNEA.IDmodel.modelParams.NB);
end
end