function [ Io ] = vech( par )
%DEVECH Inverse of the vech function
%   Vech will return a 3d symmetric matrix from its vectorization
%   Given [Ixx, Ixy, Ixz, Iyy, Iyz, Izz], returns
%   [Ixx, Ixy, Ixz; 
%    Ixy, Iyy, Iyz;
%    Ixz, Ixy, Izz]
%
    Io = zeros(3,3);
    Io(1,1) = par(1);
    Io(1,2) = par(2);
    Io(1,3) = par(3);
    Io(2,1) = par(2);
    Io(2,2) = par(4);
    Io(2,3) = par(5);
    Io(3,1) = par(3);
    Io(3,2) = par(5);
    Io(3,3) = par(6);
end

