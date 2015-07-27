function [ outputModel ] = doubleAllMasses( inputModel )
%DOUBLEALLMASSES Return a model equal to the input masses, but with double masses
outputModel = inputModel
for l = 1:length(outputModel.I)
    params = inertialParamsFromInertiaMatrix(outputModel.I{l});
    params(1) = 2*params(1)
    outputModel.I{l} = inertiaMatrixFromInertialParams(params);
end

end

