function [ outputModel ] = augmentInertialParam( inputModel, augmentFactor )
%augmentInertialParam Return a model equal to the input one, but with the
%inertial parameters augmented of the factor indicated in augmentFactor (1
%mean double, 2 means triplicate, 0.1 means an increase of 10%)
outputModel = inputModel;
for l = 1:length(outputModel.I)
    params = inertialParamsFromInertiaMatrix(outputModel.I{l});
    params = (1+augmentFactor)*params;
    outputModel.I{l} = inertiaMatrixFromInertialParams(params);
end

end

