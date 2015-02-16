function [ model ] = iCubSensLearnt( model, learn )
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


idSy_inv = []; jdSy_inv = []; dSy_inv=[];
my = 1;

if model.ny ~= length(learn.cov_ini)
   error('[ERROR] Something wrong with the learnt model, mismatch with the previous')
end

for i = 1 : model.ny
   dy = model.sizes{i,1};
   model.Sy{i,1} = learn.cov_est{i};
   
   [ii, jj, ss] = submatrixSparse(my, my, inv(model.Sy{i,1}));
   idSy_inv = [idSy_inv; ii];
   jdSy_inv = [jdSy_inv; jj];
   dSy_inv  = [dSy_inv;  ss];
   my = my + dy;
end
model.Sy_inv = sparse(idSy_inv, jdSy_inv, dSy_inv);
end
