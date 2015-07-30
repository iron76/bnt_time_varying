function [ inertialParameters ] = inertialParametersFromModel( model )
%inertialParametersFromModel Get a 10*NB vector of inertial parameters from
% a model
inertialParameters = zeros(10*model.NB,1);
for i = 1:model.NB
    inertialParameters((1+(i-1)*10):(i*10)) = inertialParamsFromInertiaMatrix(model.I{i});
end

end

