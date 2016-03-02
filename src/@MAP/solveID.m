function obj = solveID(obj)
%solveID Inverse Dynamics with sparse Newton-Euler Algorithm (SNEA)
%   This function solves the inverse dynamics problem with the sparse
%   Newton-Euler algorithm, as described in the paper "BERDY: Bayesian 
%   Estimation for Robot Dynamics. A Probabilistic Estimation of Whole-Body
%   Dynamics with Redundant Measurements." The output 'd' is structured as
%   follows:
%
%   d   = [d_1, d_2, ..., d_obj.IDstate.n]
%
%   where:
%
%   d_i = [a_i, fB_i, f_i, tau_i, fx_i, d2q_i]
%
%   and a_i is the link-i spatial acceleration, fB_i is the net spatial
%   force on the link-i, f_i is spatial wrench transmitted to link-i from
%   its parent, tau_i is torque on joint-i, fx_i is the external force on
%   link-i and d2q_i is acceleration of joint-i. The input to the algorithm
%   is in obj.IDmeas.y organized as follows:
%
%   obj.IDmeas.y = [y_1, y_2, ... , y_obj.IDsens.m]
%
%   The relationship between d and y is given by Y(q, dq) d + b_y = y where
%   the matrix Y(q, dq), is represented as a sparse matrix. Moreover, the
%   variables d should satisfy the Newton-Euler equations represented as
%   D(q,dq) d + b(q, dq) = 0, again represented as a sparse matrix. 
%
% Author: Francesco Nori
% Genova, Dec 2014

%%

NB = obj.IDmodel.modelParams.NB;
D = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);  
%D = full(D);
b_D = sparse(obj.ibs, ones(size(obj.ibs)), obj.bs, 19*NB, 1);  

% Dx = D(1:19*NB, 1:19*NB);
% Dy = D(1:19*NB, 19*NB+1:26*NB);

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
% Sw     = obj.IDmodel.modelParams.Sw.matrix;
% Sy_inv = eye(my)   ./sMeas;
%Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;
if isa(obj.IDsens.sensorsParams.Sy_inv, 'submatrixSparse')
   Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;
else
   Sy_inv = obj.IDsens.sensorsParams.Sy_inv;
end

Y = cell2mat(obj.IDsens.sensorsParams.Y);
invIndex(obj.id)= 1: length(obj.id);     %permutation vector
Y = Y(:,invIndex);

y      = obj.IDmeas.y;

S_Dinv = Sv_inv;                                %constraint equation covariance
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);  %prior covariance
S_Yinv = Sy_inv;                                %measurements equation covariance
%mu_D    = zeros(length(S_dinv), 1);
mu_d    = zeros(length(S_dinv), 1);

if isfield(obj.IDsens.sensorsParams, 'bias')
    b_Y     = obj.IDsens.sensorsParams.bias;
else
    b_Y     = zeros(size(y));
end
    

%% test with only ftx and ddq 
%  y = y(end-13:end,:);
%  Y = Y(end-13:end,:);
%  b_Y = b_Y(end-13:end,:);
%  S_Yinv = S_Yinv(end-13:end,end-13:end);
%  invIndex(obj.id)= 1: length(obj.id);
%  Y = Y(:,invIndex);
%  
%%
    
% permutations corresponding to the RNEA 
% I      = [obj.ia; obj.ifB; obj.iF(end:-1:1, 1); obj.itau];
% J      = [obj.jfx; obj.jd2q; obj.ja; obj.jfB; obj.jF(end:-1:1, 1); obj.jtau];
% Iinv(I)= 1:length(I);
% Jinv(J)= 1:length(J);
% with these definitions [Y(:,J); D(I, J)] is lowertriangular

%% mean and covariance matrix of p(d) ~ N(muBar_D, SBar_D)

SBar_D  = inv(D'*S_Dinv*D + S_dinv);
muBar_D = SBar_D*((S_dinv*mu_d) - D'*S_Dinv*b_D);

%% mean and covariance matrix of p(d|y) ~ N(mu_d|y, S_d|y)
% Note : mu_d|y --> d (as we are considering MAP estimator);
%        S_d|y  --> Sd

Sd = inv(inv(SBar_D) + Y'*S_Yinv*Y);
d  = Sd * (Y'*S_Yinv*(y-b_Y) + (inv(SBar_D) * muBar_D));

YD = [Y;D];
YD = full(YD);
YDsize = size(YD);
% d= inv([D;Y])*[-b_D;y-b_Y];

%%
obj.d  = d(obj.id,1);
obj.Sd = Sd(obj.id,obj.id);

%% ======= Nori version

% if nargin == 1
%    d      = (D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)\(Y'*S_Yinv*(y-b_Y) - D'*S_Dinv*b_D + S_dinv * mu_D);
%    obj.d  = d(obj.id,1);
%
% elseif (nargin == 2) && strcmp(computeVariance, 'variance')
%    Sd     = inv(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y);
%    d      = Sd*(Y'*S_Yinv*(y-b_Y) - D'*S_Dinv*b_D + S_dinv * mu_D);
%    obj.d  = d(obj.id,1);
%    obj.Sd = Sd(obj.id,obj.id);
% end
%    

end % solveID