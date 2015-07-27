function [ totalMass ] = getTotalMass( model )
%GETTOTALMASS Get the total mass of a spatial model
totalMass = 0.0;
for l = 1:length(model.I)
    params = inertialParamsFromInertiaMatrix(model.I{l});
    totalMass = totalMass + params(1);
end

end

