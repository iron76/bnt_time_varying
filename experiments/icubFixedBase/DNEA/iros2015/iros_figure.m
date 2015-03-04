clear all
close all
clc

res       = 0;
NB        = 10;
S_dmodel  = 1e-1;
S_ymodel  = 1e-2;


dmodel_RNEA   = autoTree(NB);
dmodel_RNEA   = autoTreeStochastic(dmodel_RNEA, S_dmodel);

ymodel_RNEA   = autoSensRNEA(dmodel_RNEA);
ymodel_RNEA   = autoSensStochastic(ymodel_RNEA, S_ymodel);

dmodel_SNEA   = dmodel_RNEA;
dmodel_SNEA.gravity = [0; -9.81; 0];

ymodel_SNEA   = autoSensSNEA(dmodel_SNEA);
ymodel_SNEA   = autoSensStochastic(ymodel_SNEA, S_ymodel);

dmodel_DNEA = dmodel_SNEA;
ymodel_DNEA = autoSensDNEA(dmodel_SNEA, ymodel_SNEA, zeros(NB,1), zeros(NB, 1));
ymodel_DNEA = autoSensStochastic(ymodel_DNEA, S_ymodel);

for i = 1 : ymodel_DNEA.ny
   lb = ymodel_DNEA.labels{i,1};
   if (length(lb) >= 7 && strcmp(lb(1:7), 'y_omega'))
      ymodel_DNEA.Sy_inv = set(ymodel_DNEA.Sy_inv, ymodel_DNEA.Sy_inv(i,i)*1e2, i, i);
      ymodel_DNEA.Sy     = set(ymodel_DNEA.Sy    , ymodel_DNEA.Sy(i,i)*1e-2   , i, i);
   end
   if (length(lb) >= 3 && strcmp(lb(1:3), 'y_q'))
      ymodel_DNEA.Sy_inv = set(ymodel_DNEA.Sy_inv, ymodel_DNEA.Sy_inv(i,i)*1e2, i, i);
      ymodel_DNEA.Sy     = set(ymodel_DNEA.Sy    , ymodel_DNEA.Sy(i,i)*1e-2   , i, i);
   end
end

dmodel_BNEA    = autoTreeStochastic(dmodel_SNEA);
ymodel_BNEA    = autoSensStochastic(ymodel_SNEA);


% res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

ymodel_RNEA = autoSensRNEA(dmodel_SNEA);
ymodel_RNEA = autoSensStochastic(ymodel_RNEA, 1e-5);

res = 0;
NB  = dmodel_SNEA.NB;

for i = 1 : 10
   q         = rand(dmodel_SNEA.NB,1);
   dq        = rand(dmodel_SNEA.NB,1);
   Sx        = diag(1e8*S_dmodel*rand(dmodel_SNEA.NB*2, 1));
   eq        = 0.8*rand(size(q));
   edq       = 1.8*rand(size(dq));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   myModel = model(dmodel_SNEA);
   mySens  = sensors(ymodel_RNEA);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   myRNEA    = RNEA(myModel, mySens);
   myRNEA    = myRNEA.setState(q,dq);
   myRNEA    = myRNEA.solveID();
   d         = myRNEA.d;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   myModel = model(dmodel_SNEA);
   mySens  = sensors(ymodel_SNEA);
   
   mySNEA    = SNEA(myModel, mySens);
   mySNEA    = mySNEA.setState(q+eq,dq+edq);
   y         = mySNEA.simY(d);
   mySNEA    = mySNEA.setY(y);
   mySNEA    = mySNEA.solveID();
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   myModel = model(dmodel_DNEA);
   mySens  = sensors(ymodel_DNEA);
   
   myDNEA    = DNEA(myModel, mySens);
   myDNEA    = myDNEA.setState(q,dq);
   y         = myDNEA.simY(d, q, dq);
   myDNEA    = myDNEA.setState(q+eq,dq+edq);
   myDNEA    = myDNEA.setY(y);
   myDNEA    = myDNEA.setStateVariance(Sx);
   d_bar     = d + rand(size(d))*100;
   myDNEA    = myDNEA.setD(d_bar);
   myDNEA    = myDNEA.solveID();
   
   
   D   = sparse(myDNEA.iDs, myDNEA.jDs, myDNEA.Ds, 19*dmodel_SNEA.NB, 26*dmodel_SNEA.NB);
   dDb = myDNEA.dDb_s.matrix;
   Ys  = myDNEA.IDsens.sensorsParams.Ys;
   DY  = [D dDb; Ys];
   
   
   if rank(full(DY)) < 26*dmodel_SNEA.NB + 2*dmodel_SNEA.NB
      disp([ 'The extended matrix [D;Y] is not full rank! Rank is: ', num2str(rank(full(DY))), ' should be ', num2str(26*dmodel_SNEA.NB + 2*dmodel_SNEA.NB)]);
      res = 1;
   else
      errSNEAd(i) = norm(d - d_bar);
      errDNEAd(i) = norm(myDNEA.d - d);
      
%       figure
%       plot(abs(mySNEA.d - d_bar), 'g')
%       hold on
%       plot(abs(myDNEA.d - d_bar))
      disp(['Diff between d.DNEA and d is ' num2str(norm(myDNEA.d - d)) ' was ' num2str(norm(d - d_bar))]);
      if norm(mySNEA.d-myDNEA.d)/length(mySNEA.d) > 0.01*max(mySNEA.d)
         disp('Result is excessively inaccurate. Test is declared failed!');
         res = 1;
      end
      
      errDNEAq(i)  = norm(myDNEA.x(   1:NB ,1) - q);
      errDNEAdq(i) = norm(myDNEA.x(NB+1:end,1) - dq);
      
      errSNEAq(i)  = norm(eq);
      errSNEAdq(i) = norm(edq);

      valDNEAq(i, :)  = abs(myDNEA.x(   1:NB ,1) - q);
      valDNEAdq(i, :) = abs(myDNEA.x(NB+1:end,1) - dq);
      
      valSNEAq(i, :)  = abs(eq);
      valSNEAdq(i, :) = abs(edq);
      
      
%       figure
%       plot(abs(myDNEA.x(1:NB,1) - q))
%       hold on
%       plot(abs(eq), 'g')
%       
%       figure
%       plot(abs(myDNEA.x(NB+1:end,1) - dq))
%       hold on
%       plot(abs(edq), 'g')
      
      disp(['Diff between x.DNEA and x was ' num2str(norm([eq; edq])) ' is ' num2str(norm(myDNEA.x-[q; dq]))]);
      if norm(myDNEA.x-[q; dq]) > norm([eq; edq])
         disp('Result is excessively inaccurate. Test is declared failed!');
         res = 1;
      end
      
   end
end

disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
disp(['Mean and variance of d-d_bar:  ' num2str(mean(errSNEAd)) '+/-' num2str(var(errSNEAd)^(1/2))]);
disp(['Mean and variance of d-d.DNEA: ' num2str(mean(errDNEAd)) '+/-' num2str(var(errDNEAd)^(1/2))]);
disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
disp(['Mean and variance of q-q_bar:  '  num2str(mean(errSNEAq))  '+/-' num2str(var(errSNEAq)^(1/2))]);
disp(['Mean and variance of q-q.DNEA: '  num2str(mean(errDNEAq))  '+/-' num2str(var(errDNEAq)^(1/2))]);
disp('\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')
disp(['Mean and variance of dq-dq_bar: ' num2str(mean(errDNEAdq)) '+/-' num2str(var(errDNEAdq)^(1/2))]);
disp(['Mean and variance of dq-dq_bar: ' num2str(mean(errSNEAdq)) '+/-' num2str(var(errSNEAdq)^(1/2))]);

figure
errorbar((1 : NB)+0.1, mean(valSNEAq,1), -2*sqrt(var(valSNEAq,1)), 2*sqrt(var(valSNEAq,1)), 'o', 'LineWidth',2, 'Color',[0.7 0.7 0.7], 'MarkerSize', 6, 'MarkerFaceColor',[0.7 0.7 0.7])
hold on;
errorbar((1 : NB)-0.1, mean(valDNEAq,1), -2*sqrt(var(valDNEAq,1)), 2*sqrt(var(valDNEAq,1)), 'o' ,'LineWidth',2, 'Color',[0 0 0], 'MarkerSize', 6, 'MarkerFaceColor',[0 0 0])
ax = gca;
axXTick = 2 : 2: NB;
for i = 2 : 2 : NB
   axXTickLabel{i/2} = strcat('q', num2str(i));
end
set(ax, 'XTick', axXTick, 'XTickLabel', axXTickLabel, 'FontSize', 24)
legend('{e}^{bar}_{q}', 'e^{map}_{q}')
xlabel('joint number')
ylabel('pos. estimation error [rad]')
grid
legend boxoff  
print -dpdf eq.pdf

figure
errorbar((1 : NB)+0.1, mean(valSNEAdq,1), -2*sqrt(var(valSNEAdq,1)), 2*sqrt(var(valSNEAdq,1)), 'o', 'LineWidth',2, 'Color',[0.7 0.7 0.7], 'MarkerSize', 6, 'MarkerFaceColor',[0.7 0.7 0.7])
hold on;
errorbar((1 : NB)-0.1, mean(valDNEAdq,1), -2*sqrt(var(valDNEAdq,1)), 2*sqrt(var(valDNEAdq,1)), 'o' ,'LineWidth',2, 'Color',[0 0 0], 'MarkerSize', 6, 'MarkerFaceColor',[0 0 0])
ax = gca;
daxXTick = 2 : 2: NB;
for i = 2 : 2 : NB
   daxXTickLabel{i/2} = strcat('dq', num2str(i));
end
set(ax, 'XTick', daxXTick, 'XTickLabel', daxXTickLabel, 'FontSize', 24)
legend('{e}^{bar}_{dq}', 'e^{map}_{dq}')
xlabel('joint number')
ylabel('vel. estimation error [rad/s]')
legend boxoff  
grid
print -dpdf edq.pdf


