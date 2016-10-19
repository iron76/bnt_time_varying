function res = testLU_ABA(dmodel)

res = 0;

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

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel_LU.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myLU    = LUABA(myModel, mySens);
myLU    = myLU.setState(q,dq);
myLU    = myLU.setY(y);
myLU    = myLU.solveID();

NB    = dmodel.NB;
tau   = zeros(NB,1);
fx    = cell(1,NB);

[P, Q, WL, WR] = factorize(myLU);

D  = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
Y  = myLU.IDsens.sensorsParams.Ys;

if norm(full(triu(P*[WL{4}*WL{3}*WL{2}*WL{1}*D*WR; Y]*Q,1))) < 1e-10
   disp('[LU_ABA] Factorisation resulted in a lower trinagular matrix')
   res = 0;
else
   disp('[LU_ABA] Error in performing the LU factorisation!!')
   res = 1;
   return
end

for i = 1 : NB
   fx{i}    = y((1:6)+(i-1)*6,1);
end
for i = 1 : NB
   tau(i,1) = y(6*NB+i);
end

[d2q_AB, pA, pa] = FDab( dmodel, q, dq, tau, fx );
[~, a_AB, fB_AB, f_AB] = ID( dmodel, q, dq, d2q_AB, fx);

d_AB = zeros(26*NB, 1);
for i = 1 : NB
   d_AB((1:26)+(i-1)*26, 1) = [a_AB{i}; fB_AB{i}; f_AB{i}; tau(i,1); fx{i}; d2q_AB(i,1)];
end

[a_LU,fB_LU,f_LU,~,~,d2q_LU] = extractDynVar(NB, myLU.d);

if res == 0 
   
   if (abs(sum(d2q_AB-d2q_LU)) > 1e-4)
      disp('[ABA] Something wrong with the d2q variable check');
      res = 1;
   end
   
   for i = 1:NB
      if (abs(sum(a_AB{1,i}-a_LU{1,i})) > 1e-4)
         disp(['[ABA] Something wrong with the a variable check ', num2str(abs(sum(a_AB{1,i}-a_LU{1,i})))]);
         res = 1;
      end
   end
   
   for i = 1:NB
      if (abs(sum(fB_AB{1,i}-fB_LU{1,i})) > 1e-4)
         disp(['[ABA] Something wrong with the fB variable check ', num2str(abs(sum(fB_AB{1,i}-fB_LU{1,i})))]);
         res = 1;        
      end
   end
   
   for i = 1:NB
      if (abs(sum(f_AB{1,i}-f_LU{1,i})) > 1e-4)
         disp(['[ABA] Something wrong with the f variable check ', num2str(abs(sum(f_AB{1,i}-f_LU{1,i})))]);
         res = 1;        
      end
   end
end

if (norm(d_AB-myLU.d) > 1e-3)
   disp(['[ABA] Something wrong with the solveID method ', num2str(norm(d_AB-myLU.d))]);
   res = 1;
end

% e_AB = residuals(myLU, d_AB)
% e_LU = residuals(myLU, myLU.d)
