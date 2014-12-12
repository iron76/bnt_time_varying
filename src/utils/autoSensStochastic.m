function [ model ] = autoSensStochastic( model )
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

sMeas   = 1e-5;

idSy_inv = []; jdSy_inv = []; dSy_inv=[];

my = 1;

for i = 1 : model.ny
   for j = 1 : model.NB
      if strcmp(model.labels{i,1}, ['a' num2str(j)])
         model.S{i,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;
      end
      if strcmp(model.labels{i,1}, ['fB' num2str(j)])
         model.S{i,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;
      end
      if strcmp(model.labels{i,1}, ['f' num2str(j)])
         model.S{i,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;
      end
      if strcmp(model.labels{i,1}, ['tau' num2str(j)])
         model.S{i,1}      = sMeas.*eye(1);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 1;
      end
      if strcmp(model.labels{i,1}, ['fx' num2str(j)])
         model.S{i,1}      = sMeas.*eye(6);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 6;
      end
      if strcmp(model.labels{i,1}, ['d2q' num2str(j)])
         model.S{i,1}      = sMeas.*eye(1);
         [ii, jj, ss] = submatrixSparse(my, my, inv(model.S{i,1}));
         idSy_inv = [idSy_inv; ii];
         jdSy_inv = [jdSy_inv; jj];
         dSy_inv  = [dSy_inv;  ss];
         my = my + 1;
      end
   end
end
model.Sy_inv = sparse(idSy_inv, jdSy_inv, dSy_inv);
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
