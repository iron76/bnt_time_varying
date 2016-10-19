clear
close all
clc

NB        = 30;
N         = 50;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

% Build the articulated chain
dmodel   = autoTree(NB, 1, 1.5, 0.9);
dmodel   = autoTreeStochastic(dmodel);
dmodel.gravity = [0; -9.81; 0];

% Build the MAP problem to solve the forward dynamics
ymodel_LU  = autoSensABA(dmodel);
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
Dc        = cell(N,1);
Yc        = cell(N,1);
Pc        = cell(N,1);
Qc        = cell(N,1);
n_fill_in = zeros(N,N);

for i = 1 : N
   q     = rand(dmodel.NB,1)*10;
   dq    = rand(dmodel.NB,1);
   myLU  = myLU.setState(q,dq);
   Dc{i} = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
   Yc{i} = myLU.IDsens.sensorsParams.Ys;
   %  Q{i} = sparse(eye(26*NB));
   [~,~,Pc{i},Qc{i}] = lu([Dc{i}; Yc{i}]);
end

for i = 1 : N
   for j = 1 : N
      % forces the LU factorisation to use diagonal pivoting
      % might result in inaccurate factorisations which
      % require the check below. Inaccurate factorisation
      % are obtained e.g. with A = [0 1; 1 0].
      [L,U,Pp] = lu(Pc{j}*[Dc{i}; Yc{i}]*Qc{j}, 0);
      if isempty(find(Pp - eye(size(Pp)), 1)) && abs(sum(sum(Pc{j}*[Dc{i}; Yc{i}]*Qc{j} - L*U))) > 1e-7
         error('unable to LU factorize')
      end
      n_fill_in(i,j) = nnz(L) + nnz(U) - nnz([Dc{i}; Yc{i}]) ;
   end
end
n_fill_in = sum(n_fill_in,1)/N;
[min_fill, index_min_fill] = min(n_fill_in);

%% 
%  Checks the computational complexity of the ABA algorithm reformulating
%  the ABA as the following factorisation: 
%
%             P_ABA*[WL*D*WR; Y]*Q_ABA 
%
%  which should result in a lower triangular matrix.
a_ABA = 205*NB - 248;
m_ABA = 224*NB - 259;   
[P_ABA, Q_ABA, WL, WR] = factorize(myLU);

D     = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*NB, 26*NB);
Y     = myLU.IDsens.sensorsParams.Ys;
b     = [sparse(myLU.ibs, ones(size(myLU.ibs)), myLU.bs, 19*NB, 1); ones(7*NB,1)];
A_ABA = tril(P_ABA*[WL{4}*WL{3}*WL{2}*WL{1}*D*WR; Y]*Q_ABA);

[a1_myABA, m1_myABA] = fw_cost(A_ABA, b);
[a2_myABA, m2_myABA] = mult_cost( WL{2}, WL{1});
[a3_myABA, m3_myABA] = mult_cost( WL{3}, WL{2}*WL{1});
[a4_myABA, m4_myABA] = mult_cost( WL{4}, WL{3}*WL{2}*WL{1});
[a5_myABA, m5_myABA] = mult_cost( WL{4}*WL{3}*WL{2}*WL{1}, D);

a_myABA = a1_myABA + a2_myABA + a3_myABA + a4_myABA + a5_myABA;
m_myABA = m1_myABA + m2_myABA + m3_myABA + m4_myABA + m5_myABA;

%% 
%  Performs a shuffle of the columns and rows so as to obtain a 
%  matrix on which we can gurantee the existence of the LU factorisation 
%  and proper defined pivots. 

B1 = D;
B2 = B1(myLU.itau, :);
B1(myLU.itau, : ) = [];
B  = [B1; B2];

A1 = [Y(:,[myLU.jfx; myLU.jtau]); B(:,[myLU.jfx; myLU.jtau])];
A2 = [Y; B];
A2(:, [myLU.jfx; myLU.jtau]) = [];
A = [A1, A2];

% Apply the LU facorisation and compute its cost.
% Pvt is the pivot matrix which is shown to coincide
% with the mass matrix 
[myL,myU,Pvt] = my_lu(A);
[a_my,m_my]   = lu_cost(A);
[a_fb ,m_fb]  = fb_cost(myL, myU, b);
a_my          = a_my + a_fb;
m_my          = m_my + m_fb;

% Apply the LU factorisation with the permutaations 
% that lead to the minimum number of fill-in. 
[minL, minU]  = lu(Pc{index_min_fill}*[D; Y]*Qc{index_min_fill}, 0);
[a_min,m_min] = lu_cost(A);
[a_fb ,m_fb]  = fb_cost(minL, minU, b);
a_min         = a_min + a_fb;
m_min         = m_min + m_fb;

X = cell(NB,NB);
for i = 1 : NB
   X{i,i} = eye(6);
   for j = i+1 : NB
      X{j,i} = myLU.Xup{j} * X{j-1,i};
   end
end

% Check that R = Pvt, being:
%
% For j >= i we have: Rij = sum_{k = j}^NB S_i' X_{k,i}' I_k X_{k,j} S_j 
%
R = zeros(NB,NB);
for i = 1 : NB
   for j = i : NB
      for k = j : NB
         R(i,j) = R(i,j) + myLU.IDmodel.S{i}' * X{k,i}' * myLU.IDmodel.modelParams.I{k} *  X{k,j} * myLU.IDmodel.S{j};
      end
      R(j,i) = R(i,j);
   end
end

if norm(full(Pvt{end-NB+1}) - R) > 1e-10
   error('error in computing the partial LU factor');
end

% Check that the mass matrix H equals Pvt=R
%
[H,C] = HandC( dmodel, q, dq);
if norm(H - R) > 1e-10
   error('error in computing the partial LU factor');
end

