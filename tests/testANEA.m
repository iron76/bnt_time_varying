function res = testANEA(dmodel, ymodel, dmodel_ANEA, ymodel_ANEA)

res = 0;

% clear all
% close all
% clc
% 
% res       = 0;
% NB        = 30;
% S_dmodel  = 1e-2;
% S_ymodel  = 1e-4;
% 
% dmodel         = autoTree(NB);
% dmodel    = autoTreeStochastic(dmodel, S_dmodel);
% dmodel.gravity = [0; -9.81; 0];
% 
% ymodel   = autoSensRNEA(dmodel);
% ymodel   = autoSensStochastic(ymodel, S_ymodel);
% 

NB        = dmodel.NB;
q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y_RNEA    = rand(ymodel.m,1);
y_ANEA    = y_RNEA(7:7:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(model(dmodel), sensors(ymodel));
myRNEA    = myRNEA.setState(q,dq);
myRNEA    = myRNEA.setY(y_RNEA);
myRNEA    = myRNEA.solveID();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myANEA    = ANEA(model(dmodel_ANEA), sensors(ymodel_ANEA));
myANEA    = myANEA.setState(q,dq);
myANEA    = myANEA.setY(y_ANEA);
myANEA    = myANEA.solveID();

if (sum(q-myANEA.IDstate.q))
   disp('Something wrong with the setQ method');
   res = 1;
end

if (sum(dq-myANEA.IDstate.dq))
   disp('Something wrong with the setDq method');
   res = 1;
end

if (sum(y_ANEA-myANEA.IDmeas.y))
   disp('Something wrong with the setY method');
   res = 1;
end

[a_RNEA,~,~,~,~,d2q_RNEA] = extractDynVar(NB, myRNEA.d);

for i = 1 : NB
   a_ANEA   = myANEA.d((i-1)*7+1 : (i-1)*7+6, 1);
   d2q_ANEA = myANEA.d((i-1)*7+7 : (i-1)*7+7, 1);
   % disp(['Diff a between RNEA and ANEA is ' num2str(norm(a_ANEA-a_RNEA{i}))]);
   if norm(a_ANEA-a_RNEA{i}) > 0.01
      disp(['Result a_ANEA-a_RNEA is excessively inaccurate. Test is declared failed! Error is: ' num2str(norm(a_ANEA-a_RNEA{i}))]);
      res = 1;
      
   end
   % disp(['Diff d2q between RNEA and ANEA is ' num2str(norm(d2q_ANEA-d2q_RNEA(i)))]);
   if norm(d2q_ANEA-d2q_RNEA(i)) > 0.01
      disp(['Result a_ANEA-a_RNEA is excessively inaccurate. Test is declared failed! Error is: ' num2str(norm(d2q_ANEA-d2q_RNEA(i)))]);
      res = 1;      
   end
end


num_of_tests  = 1;
ymodel_DANEA = autoSensDANEA(dmodel_ANEA, ymodel_ANEA, zeros(NB,1), zeros(NB, 1));
myANEA = DANEA(model(dmodel_ANEA), sensors(ymodel_DANEA));


for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 7,1);
   
   myANEA = myANEA.setState(q,dq);
   
   f  = @(x) computeY(myANEA, d, x);
   dy = deriv(f, [q; dq]);
   
   % In current formulation 
   %
   % y = Yd d + Yx x + b_Y(x) = [Yd(x) Yx][d ; x] + b_Y(x)
   %
   % and therefore the derivative is:
   %
   % dy/dx = Yx + db_Y/dx
   
   Yx  = myANEA.IDsens.sensorsParams.Ys(:,(end-2*dmodel.NB+1):end);
   dby = myANEA.dby_s.matrix + Yx;
   if norm( dby - dy) > 1e-6
      disp(['dy numerical derivative is quite different: ' num2str(norm(dby - dy))])
      imagesc([dby, dy])
      colorbar
      res = 1;
   end
end


for i = 1 : num_of_tests
   q  = rand(dmodel.NB     ,1);
   dq = rand(dmodel.NB     ,1);
   d  = rand(dmodel.NB * 7,1);
   Sx = diag(rand(dmodel.NB*2, 1));
   f  = @(x) computeDb(myANEA , d , x);
   dD = deriv(f, [q; dq]);
   
   myANEA = myANEA.setState(q,dq);
   myANEA = myANEA.setD(d);
   myANEA = myANEA.setDprior(d);
   myANEA = myANEA.setXprior([q; dq]);
   myANEA = myANEA.setXvariance(Sx);
   if norm(myANEA.dDb_s.matrix - dD) > 1e-3
      disp(['dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb_s.matrix - dD))])
      res = 1;
   end

   if norm(myANEA.dDb.matrix - dD) > 1e-3
      disp(['dD numerical derivative is quite different: ' num2str(norm(myDNEA.dDb.matrix - dD))])
      res = 1;
   end
   
   if norm(full(myANEA.Sx_inv.matrix) - inv(Sx)) ~= 0
      disp('Sx_inv was not properly set!')
      res = 1;
   end
   
   if norm(full(myANEA.Sx.matrix) - Sx) ~= 0
      disp('Sx was not properly set!')
      res = 1;
   end

end


