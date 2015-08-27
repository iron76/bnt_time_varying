function [ D ] = computeD( myNEA , q)

[~,n] = size(q);
if n ~= 1
   error('In computing D(x), x should be a column vector')
end

dq     = zeros(size(q));
myNEA  = myNEA.setState(q, dq);
D      = sparse(myNEA.iDs, myNEA.jDs, myNEA.Ds, 19*myNEA.IDmodel.modelParams.NB, 26*myNEA.IDmodel.modelParams.NB);
end