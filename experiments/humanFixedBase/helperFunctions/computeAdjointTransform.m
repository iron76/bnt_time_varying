function [ X_B_A ] = computeAdjointTransform(R_B_A, r_A_from0toP)
%COMPUTEADJOINTRANSFORM return B_X_A as in Featherstone notation 
% [Featherstone(2008), Rigid Body Dynamics Algorithms (2008)].
% Considering:
%  - a reference frame A with origin in 0;
%  - a reference frame B with origin in P;
%
% Inputs:   
%  - B_R_A          :rotation matrix from A to B;
%  - A_r_from0toP   :3x1 position vector from O to P expressed in A;

X_B_A =   [          R_B_A                  zeros(3) ; 
           -R_B_A*skew(r_A_from0toP)         R_B_A   ];
        
end

