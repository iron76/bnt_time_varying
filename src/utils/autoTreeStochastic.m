function [ model ] = autoTreeStochastic( model )
%AUTOTREESTOCHASTIC Add stochastic component to a dynamic model
%   This function takes a structure containing an articulated rigid body
%   model (similar to the one created by autoSensSNEA) and adds to the
%   structure some fields that are used to represent the variance of
%   Newton-Euler equations describing the model itself. This function is
%   used to implement some of the ideas described in "BERDY: Bayesian
%   Estimation for Robot Dynamics. A Probabilistic Estimation of Whole-Body
%   Dynamics with Redundant Measurements." The Newton-Euler equations are
%   represented as follows:
%
%                  D(q,dq) d + b(q, dq) = 0               (1)
%
%   D has dimension (19NB, 26NB) being NB the number of rigid bodies in the
%   articulated chain. The 19NB rows are subvided in NB groups of 19 rows
%   and the 19 rows subdivided in groups of 6, 6, 6, and 1. Each of these
%   sub-sub-groups has an associated variance, denoted Sm.a, Sm.fB, Sm.f,
%   and Sm.tau respectively. The global variance for equation (1) is
%   represented in a sparse matrix Sv and its inverse stored in sparse
%   model.Sv_inv. Additional information is stored as a prior on the vector
%   d which has the following structure:
%
%              d_j = [a_j, fB_j, f_j, tau_j, fx_j, d2q_j]
%
%                d = [d_1, ... , d_NB]
%
%  A variance for the prior of fx_j, d2q_j is also represented in the matrices
%  Su.fx and Su.d2q, and the overall covariance matrix Sw represented as a
%  sparse matrix in its inverse model.Sw_inv

sModel  = 1e-5;
sUknown = 1e3;

idSw_inv = []; jdSw_inv = []; dSw_inv=[];
idSv_inv = []; jdSv_inv = []; dSv_inv=[];

for i = 1 : model.NB
   model.Sm.a{i}   = sModel.*generateSPDmatrix(6);
   model.Sm.fB{i}  = sModel.*generateSPDmatrix(6);
   model.Sm.f{i}   = sModel.*generateSPDmatrix(6);
   model.Sm.tau{i} = sModel.*generateSPDmatrix(1);
   
   [ii, jj, ss] = submatrixSparse(19*i-18, 19*i-18, ...
      [inv(model.Sm.a{i})  zeros(6,6)    zeros(6,6)   zeros(6,1); ...
      zeros(6,6)     inv(model.Sm.fB{i}) zeros(6,6)   zeros(6,1); ...
      zeros(6,6)     zeros(6,6)    inv(model.Sm.f{i}) zeros(6,1); ...
      zeros(1,6)     zeros(1,6)    zeros(1,6)   inv(model.Sm.tau{i})]);
   idSv_inv = [idSv_inv; ii];
   jdSv_inv = [jdSv_inv; jj];
   dSv_inv  = [dSv_inv;  ss];
   
   model.Su.fx{i}  = sUknown.*generateSPDmatrix(6);
   model.Su.d2q{i} = sUknown.*generateSPDmatrix(1);
   
   [ii, jj, ss] = submatrixSparse(7*i-6, 7*i-6, [inv(model.Su.fx{i}) zeros(6,1); zeros(1,6) inv(model.Su.d2q{i})]);
   idSw_inv = [idSw_inv; ii];
   jdSw_inv = [jdSw_inv; jj];
   dSw_inv  = [dSw_inv;  ss];
   
   
end
model.Sv_inv = sparse(idSv_inv, jdSv_inv, dSv_inv);
model.Sw_inv = sparse(idSw_inv, jdSw_inv, dSw_inv);
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

