function [ model ] = autoSensStochastic( model , sMeas)
%AUTOSENSSTOCHASTIC Add stochastic component to a sensor distribution
%   This function takes a structure containing sensor distribution (similar
%   to the one created by autoSensSNEA) and adds to the structure some
%   fields that are used to represent the variance of the measurements. The
%   sensor model is the following:
%
%                 Y(q,dq) d = y
%
%   and the variance is associated with the confidence on the measurement
%   equation.

if nargin == 1
   sMeas     = 1;
   generateS = @(n)eye(n);
else
   generateS = @(n)generateSPDmatrix(n);
end

iSy_s = cell2mat(model.sizes);
jSy_s = cell2mat(model.sizes);
model.Sy_inv = submatrixSparse(iSy_s, jSy_s, (1:length(iSy_s))', (1:length(jSy_s))');
model.Sy     = submatrixSparse(iSy_s, jSy_s, (1:length(iSy_s))', (1:length(jSy_s))');

for i = 1 : model.ny
   dy = model.sizes{i,1};
   % model.Sy{i,1} = sMeas.*generateSPDmatrix(dy);
   % S = inv(model.Sy{i,1});
   S = sMeas.*generateS(dy);
   model.Sy_inv = set(model.Sy_inv, inv(S), i, i);
   model.Sy     = set(model.Sy    ,     S , i, i);
end

end

function A = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix

A = rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix using either
A = A+A';
% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
A = A + n*eye(n);

end
