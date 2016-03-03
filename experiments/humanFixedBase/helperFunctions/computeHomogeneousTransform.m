function [ H_A_B ] = computeHomogeneousTransform(R_A_B, r_A_fromAtoB)
%COMPUTEHOMOGENEOUSTRANSFORM return A_H_B.
% Considering:
%  - a reference frame A with origin in 0_A;
%  - a reference frame B with origin in 0_B;
%  - a point P originally expressed in frame B that we want to express in
%    frame A
%
% Inputs:   
%  - A_R_B          :rotation matrix from B to A;
%  - A_r_fromAtoB   :3x1 position vector from A to B expressed in A;

H_A_B =   [  R_A_B       r_A_fromAtoB ; 
               0                1     ];
        
end

