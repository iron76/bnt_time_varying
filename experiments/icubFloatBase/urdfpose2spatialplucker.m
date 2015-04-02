function [ X ] = urdfpose2spatialplucker( rpy, xyz )
%URDFPOSE2SPATIALPLUCKER Transform a urdf pose 
%   into a plucker transformation 
     R = rotx(rpy(1))*roty(rpy(2))*rotz(rpy(3));
     R = R(1:3,1:3);
     X = pluckerFromSE3Transform(R,[0.0 0.0 0.0]);
end

