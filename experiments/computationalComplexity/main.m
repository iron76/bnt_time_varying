clear
close all
clc

NB        = 30;
N         = 50;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

dmodel   = autoTree(NB, 1, 1.5, 0.9);
dmodel   = autoTreeStochastic(dmodel);
dmodel.gravity = [0; -9.81; 0];

ymodel_LU  = autoSensABA(dmodel);
ymodel_LU  = autoSensStochastic( ymodel_LU );
mySens  = sensors( ymodel_LU );

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
% % set Sy_inv to zero
% i = mySens.sensorsParams.Sy_inv.i;
% j = mySens.sensorsParams.Sy_inv.j;
% for k = 1:length(i)
%    mySens.sensorsParams.Sy_inv = set(mySens.sensorsParams.Sy_inv, eye(size(mySens.sensorsParams.Sy_inv(i(k),j(k)))), i(k), j(k));
% end

myModel = model(dmodel_LU);
myLU    = LUABA(myModel, mySens);
D       = cell(N,1);
Y       = cell(N,1);
P       = cell(N,1);
Q       = cell(N,1);
nzLU    = zeros(N,N);
nzA     = zeros(N,N);

for i = 1 : N
   q    = rand(dmodel.NB,1)*10;
   dq   = rand(dmodel.NB,1);
   myLU = myLU.setState(q,dq);
   D{i} = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
   Y{i} = myLU.IDsens.sensorsParams.Ys;
   %  Q{i} = sparse(eye(26*NB));
   [L,U,P{i},Q{i}] = lu([D{i}; Y{i}]);
end

for i = 1 : N
   for j = 1 : N
      % forces the LU factorisation to use diagonal pivoting
      % might result in inaccurate factorisations which
      % require the check below. Inaccurate factorisation
      % are obtained e.g. with A = [0 1; 1 0].
      [L,U,Pp] = lu(P{j}*[D{i}; Y{i}]*Q{j}, 0);
      if isempty(find(Pp - eye(size(Pp)), 1)) && abs(sum(sum(P{j}*[D{i}; Y{i}]*Q{j} - L*U))) < 1e-8
         nzLU(i,j) = nnz(L) + nnz(U);
         nzA(i,j) = nnz(P{j}*[D{i}; Y{i}]*Q{j});
      else
         error('unable to LU factorize')
      end
   end
end

% plot(sum(nzLU,1))

[P, Q, WL, WR] = factorize(myLU);
for i = 1 : length(myLU.jd2q)
   jd2q_a(1+7*(i-1):7*i,1) = [myLU.jd2q(i); myLU.ja(1+6*(i-1):6*i)];
end
for i = 1 : length(myLU.jtau)
   itau_a(1+7*(i-1):7*i,1) = [myLU.itau(i); myLU.ia(1+6*(i-1):6*i)];
end
Dt = D{N} * WR;
Yt = Y{N} * WR;
pt = [ myLU.ifB(end:-1:1, 1);  itau_a(end:-1:1);  myLU.iF;];
qt = [ myLU.jfB(end:-1:1, 1);  jd2q_a(end:-1:1);  myLU.jF; myLU.jfx; myLU.jtau;];
A  = [Dt(pt,qt); Yt(:,qt)];

[L,U] = my_lu(A);
pivot = diag(U);
for i = 1 : NB
   j = NB-i+1;
   if pivot(6*(NB+1)+1+7*(i-1)) - myLU.IDmodel.S{j}'*myLU.IA{j}*myLU.IDmodel.S{j} ~=0
      error('error in checking the pivot');
   end
end

B1 = D{N};
B2 = B1(myLU.itau, :);
B1(myLU.itau, : ) = [];
B  = [B1; B2];

A1 = [Y{N}(:,[myLU.jfx; myLU.jtau]); B(:,[myLU.jfx; myLU.jtau])];
A2 = [Y{N}; B];
A2(:, [myLU.jfx; myLU.jtau]) = [];
A = [A1, A2];

[L,U,P] = my_lu(A);
pivot = diag(U);

X = cell(NB,NB);
for i = 1 : NB
   X{i,i} = eye(6);
   for j = i+1 : NB
      X{j,i} = myLU.Xup{j} * X{j-1,i};
   end
end

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

if norm(full(P{end-NB+1}) - R) > 1e-10
   error('error in computing the partial LU factor');
end

[H,C] = HandC( dmodel, q, dq);
if norm(H - R) > 1e-10
   error('error in computing the partial LU factor');
end


showmotion( dmodel, [1 2], repmat(q,1,2));

% [L, U, P, Q] = lu(A);
% for i = 1 : N
%    I = find(P(:,i));
%    J = find(Q(i,:));
%    disp(['[' num2str(i) ',' num2str(i) '] -> ' '[' num2str(I) ',' num2str(J) ']' ])
% end
% subplot(122)
% imagesc(P*A*Q);


% for i = 1 : N
%    dP = P{1} - P{i};
%    [I,J] = find(dP);
%    for k = 1 : length(J)
%       disp([myLU.iLabels(J(k)) ' -> ' num2str(I(k))])
%    end
% end
% 
% for i = 1 : N
%    dQ = Q{1} - Q{i};
%    [I,J] = find(dQ);
%    for k = 1 : length(I)
%       disp([myLU.jLabels(I(k)) ' -> ' num2str(J(k))])
%    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% z = rand(6,1);
% y = rand(6,1);
% (z'*myLU.IDmodel.modelParams.I{2}*z)*(y'*myLU.IDmodel.modelParams.I{2}*y) - (y'*myLU.IDmodel.modelParams.I{2}*z)^2