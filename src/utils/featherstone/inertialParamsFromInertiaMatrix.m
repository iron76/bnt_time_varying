function [ inertialParams ] = inertialParamsFromInertiaMatrix( inertiaMatrix )
%inertialParamsFromInertiaMatrix Convert an inertial matrix to a vector of
% 10 parameters
% Serialization assumed : m mcx mcy mcz Ixx Ixy Ixz Iyy Iyz Izz
inertialParams = zeros(10,1);

mass = inertiaMatrix(6,6);
% apparently skew is both skew and its inverse function... a bit confusing. 
mcom = skew(inertiaMatrix(1:3,4:6));
inertiaAtFrameOrigin = inertiaMatrix(1:3,1:3);

inertialParams(1) = mass;
inertialParams(2:4) = mcom;
inertialParams(5:10) = vech(inertiaAtFrameOrigin);

end

