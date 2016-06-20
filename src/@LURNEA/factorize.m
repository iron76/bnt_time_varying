function [P, Q] = factorize(obj)

% Author: Francesco Nori
% Genova, June 2016
%
% Given an inverse-dynamics problem formulated as [D; Y]*d + [bD; -y] = 0
% with d organized as follows:
%
%   d   = [dx_1, dx_2, ..., dx_obj.IDstate.n, dy_1, dy_2, ..., dy_obj.IDstate.n]
%
%   where:
%
%   dx_i = [a_i, fB_i, f_i, tau_i]
%
%   dy_i = [fx_i, d2q_i]
%
% computes permutation matrices P and Q such that P*[D;Y]*Q = L with L
% lower triangular. 

NB = obj.IDmodel.modelParams.NB;

p  = [obj.ia; obj.ifB; obj.iF(end:-1:1, 1); obj.itau];
q  = [obj.jfx; obj.jd2q; obj.ja; obj.jfB; obj.jF(end:-1:1, 1); obj.jtau];
D  = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
Y  = obj.IDsens.sensorsParams.Ys;
I  = eye(size([D;Y]));
Q  = I(:,q);

ID = eye(size(D*D'));
IY = eye(size(Y*Y'));
P  = blkdiag(IY, ID(p,:))*[0*Y*D', IY; ID, 0*D*Y'];

% The resulting P*[D; Y]*Q is lower triangular

return




