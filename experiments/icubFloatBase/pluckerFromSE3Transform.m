function [ X ] = pluckerFromSE3Transform( R, xyz )
%UNTITLED Transform {}^B T_A to {}^B X_A
X = [R, zeros(3,3); skew(xyz)*R, R];
 
end

