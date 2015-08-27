function [ inertiaMatrix ] = inertiaMatrixFromInertialParams( inertialParams )
%inertiaMatrixFromInertialParams Convert a vector of 10 inertial params to
%an inertial Matrix
% Serialization assumed : m mcx mcy mcz Ixx Ixy Ixz Iyy Iyz Izz
inertiaMatrix = zeros(6,6);

mass = inertialParams(1);
mcom = inertialParams(2:4);
inertiaAtFrameOrigin = devech(inertialParams(5:10));

inertiaMatrix(1:3,1:3) = inertiaAtFrameOrigin;
inertiaMatrix(1:3,4:6) = skew(mcom);
inertiaMatrix(4:6,1:3) = -skew(mcom);
inertiaMatrix(4:6,4:6) = mass*eye(3);

end

