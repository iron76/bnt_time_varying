function res = testLU_RNEA(dmodel)

res = 0;

ymodel  = autoSensRNEA(dmodel);
ymodel  = autoSensStochastic( ymodel );
mySens  = sensors( ymodel );

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
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myLU    = LURNEA(myModel, mySens);
myLU    = myLU.setState(q,dq);
myLU    = myLU.setY(y);
myLU    = myLU.solveID();

NB    = dmodel.NB;
d2q   = zeros(NB,1);
fx    = cell(1,NB);

[P, Q] = factorize(myLU);

D  = sparse(myLU.iDs, myLU.jDs, myLU.Ds, 19*dmodel.NB, 26*dmodel.NB);
Y  = myLU.IDsens.sensorsParams.Ys;

if isempty(find(tril(P*[D; Y]*Q)-P*[D; Y]*Q, 1))
   disp('[LU_RNEA] Factorisation resulted in a lower trinagular matrix')
   res = 0;
else
   disp('[LU_RNEA] Error in performing the LU factorisation!!')
   res = 1;
   return
end

for i = 1 : NB
   fx{i}    = y((1:6)+(i-1)*6,1);
end
for i = 1 : NB
   d2q(i,1) = y(6*NB+i); 
end

[tau_ID, a_ID, fB_ID, f_ID] = ID( dmodel, q, dq, d2q, fx);

d_ID = zeros(26*NB, 1);
for i = 1 : NB
   d_ID((1:26)+(i-1)*26, 1) = [a_ID{i}; fB_ID{i}; f_ID{i}; tau_ID(i,1); fx{i}; d2q(i,1)];
end

[a_LU,fB_LU,f_LU,tau_LU,~,~] = extractDynVar(NB, myLU.d);

if res==0
   if (abs(sum(tau_LU-tau_ID)) > 1e-4)
      disp(['[LU_RNEA] Something wrong with the tau variable check ', num2str(abs(sum(tau_LU-tau_ID)))]);
      res = 1;
   end

   for i = 1:NB
      if (abs(sum(f_LU{1,i}-f_ID{1,i})) > 1e-4)
         disp(['[LU_RNEA] Something wrong with the f variable check ', num2str(abs(sum(f_LU{1,i}-f_ID{1,i})))]);
         res = 1;
      end
   end
      
   for i = 1:NB
      if (abs(sum(a_LU{1,i}-a_ID{1,i})) > 1e-4)
         disp(['[LU_RNEA] Something wrong with the a variable check ', num2str(abs(sum(a_LU{1,i}-a_ID{1,i})))]);
         res = 1;
      end
   end
   
   for i = 1:NB
      if (abs(sum(fB_LU{1,i}-fB_ID{1,i})) > 1e-4)
         disp(['[LU_RNEA] Something wrong with the fB variable check ', num2str(abs(sum(fB_LU{1,i}-fB_ID{1,i})))]);
         res = 1;
      end
   end      
end

if (norm(d_ID-myLU.d) > 1e-3)
   disp(['[LU_RNEA] Something wrong with the solveID method ', num2str(norm(d_ID-myLU.d))]);
   res = 1;
end

% e_ID = residuals(myLU, d_ID)
% e_LU = residuals(myLU, myLU.d)