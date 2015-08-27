function [ regressor ] = inertia3DRegressor( omega )
%INERTIA3DREGRESSOR Returns the 6x10 inertia regressor 
%   Returns the 3x6 inertia regressor such that:
%   inertia3DRegressor(omega)*vech(Io) ==
%   Io*omega
    wx = omega(1);
    wy = omega(2);
    wz = omega(3);
    regressor = zeros(3,6);
    
    % first line
    regressor(1,1) = wx;
    regressor(1,2) = wy;
    regressor(1,3) = wz;
    
    % second line 
    regressor(2,2) = wx;
    regressor(2,4) = wy;
    regressor(2,5) = wz;
    
    % third line 
    regressor(3,3) = wx;
    regressor(3,5) = wy;
    regressor(3,6) = wz;

end

