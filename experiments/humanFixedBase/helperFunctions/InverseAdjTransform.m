function [B] = InverseAdjTransform(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

block1 = A(1:3,1:3)';
block2 = A(1:3,4:6)';
block3 = A(4:6,1:3)';
block4 = A(4:6,4:6)';

B = [block1  block2;
     block3  block4];
 
end
