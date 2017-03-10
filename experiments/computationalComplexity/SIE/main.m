clear
close all
clc

NB        = 12;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

% Build the articulated chain
dmodel   = autoTree(NB, 1, 1.5, 1);
dmodel   = autoTreeStochastic(dmodel);
dmodel.gravity = [0; -9.81; 0];

% Build the MAP problem to solve the forward dynamics
ymodel_LU  = autoSensSIEfriction(dmodel);
ymodel_LU  = autoSensStochastic( ymodel_LU );
mySens     = sensors( ymodel_LU );
dmodel_LU  = autoTreeStochastic(dmodel);
% set Sw_inv to zero
i = dmodel_LU.Sw_inv.i;
j = dmodel_LU.Sw_inv.j;
for k = 1:length(i)
   dmodel_LU.Sw_inv = set(dmodel_LU.Sw_inv, zeros(size(dmodel_LU.Sw_inv(i(k),j(k)))), i(k), j(k));
end
% set Sv_inv to zero
i = dmodel_LU.Sv_inv.i;
j = dmodel_LU.Sv_inv.j;
for k = 1:length(i)
   dmodel_LU.Sv_inv = set(dmodel_LU.Sv_inv, eye(size(dmodel_LU.Sv_inv(i(k),j(k)))), i(k), j(k));
end

%% 
%  Numerically check if there exist P and Q in the LU factorisation:
%
%             P*[D(q);Y(q)]*Q = L*U 
%
%  which do not depend on thespecific q.

myModel   = model(dmodel_LU);
myLU      = LUABA(myModel, mySens);

q     = rand(dmodel.NB,1)*10;
dq    = rand(dmodel.NB,1);
myLU  = myLU.setState(q,dq);
D     = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
Y     = myLU.IDsens.sensorsParams.Ys;
b     = [sparse(myLU.ibs, ones(size(myLU.ibs)), myLU.bs, 19*NB, 1); ones(7*NB,1)];

% Apply the LU factorisation with the permutaations
% that lead to the minimum number of fill-in.
[~,~,P,Q] = lu([D; Y]);

[minL, minU]  = lu(P*[D; Y]*Q, 0);
[a1_min,m1_min] = lu_cost(P*[D; Y]*Q);
[a1_fb ,m1_fb]  = fb_cost(minL, minU, b);
a_PLUQABA(NB)     = a1_min + a1_fb;
m_PLUQABA(NB)     = m1_min + m1_fb;

showmotion( dmodel, [1 2], repmat(q,1,2));

