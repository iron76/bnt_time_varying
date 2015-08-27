function [ regressor ] = netWrenchRegressor( a , v )
%NETWRENCHGREGRESSOR Returns the 6x10 netWrenchRegressor
%   Returns the 6x10 matrix such that:
%   netWrenchRegressor(a,v)*inertialParamsFromInertiaMatrix(I) ==
%   I*a + crf(v)*I*v
regressor = inertiaRegressor(a) + crf(v)*inertiaRegressor(v);

end

