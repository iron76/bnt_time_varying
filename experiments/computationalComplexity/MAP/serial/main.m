clear
close all
clc

NB        = 30;
N         = 50;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

% Build the articulated chain
dmodel   = autoTree(NB, 1, 1, 1);
dmodel   = autoTreeStochastic(dmodel);
dmodel.gravity = [0; -9.81; 0];

% Build the MAP problem to solve the forward dynamics
ymodel_LU  = autoSensABA(dmodel);
ymodel     = addSens(ymodel_LU, 'tau2');
ymodel     = autoSensStochastic( ymodel );
mySens     = sensors( ymodel );
dmodel     = autoTreeStochastic(dmodel);
% set Sw_inv to zero
i = dmodel.Sw_inv.i;
j = dmodel.Sw_inv.j;
for k = 1:length(i)
   dmodel.Sw_inv = set(dmodel.Sw_inv, zeros(size(dmodel.Sw_inv(i(k),j(k)))), i(k), j(k));
end
% set Sv_inv to zero
i = dmodel.Sv_inv.i;
j = dmodel.Sv_inv.j;
for k = 1:length(i)
   dmodel.Sv_inv = set(dmodel.Sv_inv, eye(size(dmodel.Sv_inv(i(k),j(k)))), i(k), j(k));
end

%% 
%  Numerically check if there exist P and Q in the LU factorisation:
%
%             P*[D(q);Y(q)]*Q = L*U 
%
%  which do not depend on thespecific q.

myModel = model(dmodel);
myMAP   = MAP(myModel, mySens);
q       = rand(dmodel.NB,1)*10;
dq      = rand(dmodel.NB,1);
y       = rand(ymodel.m,1);
myMAP   = myMAP.setState(q,dq);
myMAP   = myMAP.setY(y);
myMAP   = myMAP.solveID();


D = sparse(myMAP.iDs, myMAP.jDs, myMAP.Ds, 19*NB, 26*NB);
Y = myMAP.IDsens.sensorsParams.Ys;
b = [-sparse(myMAP.ibs, ones(size(myMAP.ibs)), myMAP.bs, 19*NB, 1); y];
d = (D'*D + Y'*Y)\([D' Y']*b);
d = d(myMAP.id,1);

if norm(d-myMAP.d) > 1e-12
   error(['The MAP computations are not consistent with the least-squares solution. Error is: ', num2str(norm(d-myMAP.d))])
end

%% 
%  Checks the computational complexity of the ABA algorithm reformulating
%  the ABA as the following factorisation: 
%
%             P_ABA*[WL*D*WR; Y]*Q_ABA 
%
%  which should result in a lower triangular matrix.
a_ABA = 205*NB - 248;
m_ABA = 224*NB - 259;   


% Apply the LU facorisation and compute its cost.
% Pvt is the pivot matrix which is shown to coincide
% with the mass matrix 
[U,p,S] = chol(D'*D + Y'*Y);
if p ~= 0
   error('Unable to perform the Cholesky factorisation');
end

L = my_chol(S'*(D'*D + Y'*Y)*S);
[ ach, mch ] = chol_cost( S'*(D'*D + Y'*Y)*S );
[ afb, mfb ] = fb_cost( L, L', [D' Y']*b );

a_map = ach + afb;
m_map = mch + mfb;


showmotion( dmodel, [1 2], repmat(q,1,2));


