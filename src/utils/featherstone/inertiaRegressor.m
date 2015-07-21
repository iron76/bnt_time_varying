function [ regressor ] = inertiaRegressor( v )
%INERTIAREGRESSOR Returns the 6x10 inertia regressor 
%   Returns the 6x10 inertia regressor such that:
%   inertiaRegressor(v)*inertialParamsFromInertiaMatrix(I) ==
%   I*v
    vLin = v(4:6);
    vAng = v(1:3);
    regressor = zeros(6,10);
    
    regressor(1:3,1)   = zeros(3,1);
    regressor(1:3,2:4) = -skew(vLin);
    regressor(1:3,5:10) = inertia3DRegressor(vAng);
    
    regressor(4:6,1)   = vLin;
    regressor(4:6,2:4) = skew(vAng);
    regressor(4:6,5:10) = zeros(3,6);

end

