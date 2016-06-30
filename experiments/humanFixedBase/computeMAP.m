function [mu_dgiveny, Sigma_dgiveny, mapObj ] = computeMAP( mapObj,q,dq,Ymatrix, y,b_y, Sigma_ygivend )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


% 1. check initial parameters;

% 2. function;
% setSigma_ygivend in MAP

mapObj = mapObj.setState(q, dq);
mapObj = mapObj.setY(y);
mapObj = mapObj.setYmatrix(Ymatrix);
mapObj = mapObj.setBias(b_y);
mapObj = mapObj.setMeasCovariance(Sigma_ygivend);

mapObj = mapObj.solveID();

% 3. formatting output
mu_dgiveny = mapObj.d;
Sigma_dgiveny = mapObj.Sd;

end