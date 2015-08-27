function [ vechI ] = vech( I )
%VECH Vectorize a symmetric matrix
%   Vech will return the vectorization of a 3d symmetric matrix 
%   Returns [Ixx, Ixy, Ixz, Iyy, Iyz, Izz]
    vechI = zeros(6,1);
    vechI(1) = I(1,1);
    vechI(2) = I(1,2);
    vechI(3) = I(1,3);
    vechI(4) = I(2,2);
    vechI(5) = I(2,3);
    vechI(6) = I(3,3);
end

