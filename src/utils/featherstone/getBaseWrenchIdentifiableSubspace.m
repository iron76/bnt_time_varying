function [ identifiableBasis , notIdentifiableBasis ] = getBaseWrenchIdentifiableSubspace( model )
%getBaseWrenchIdentifiableSubspace Get the B matrix such that 
%  baseParam = B'*inertialParams, i.e.  the matrix whose columns 
% are a base of the identifiable parameters

% implement gautier algorithm
% 

nrOfSamples = 100;
inertialParamSize = model.NB*10;

A = zeros(inertialParamSize,inertialParamSize);

for i = 1:nrOfSamples 
    q = rand(model.NB,1);
    dq = rand(model.NB,1);
    ddq = rand(model.NB,1);
    
    fBaseRegr = regressorID(model, q, dq, ddq);
    
    A = A + fBaseRegr'*fBaseRegr;
end

[U,S,V] = svd(A);

baseParamSize = rank(A);

identifiableBasis    = V(:, 1:baseParamSize);
notIdentifiableBasis = V(:, (baseParamSize+1):inertialParamSize);

end