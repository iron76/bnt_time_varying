clear all
close all
clc

NB      = 30;
n       = 10;

dmodel  = autoTree(NB);
ymodel  = autoSensSNEA(dmodel);

dmodel  = autoTreeStochastic(dmodel);
ymodel  = autoSensStochastic(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myBNEA  = BNEA(myModel, mySens);

bnet    = cell(1,n);
engine  = cell(1,n);
sample  = cell(NB*6+myBNEA.IDsens.sensorsParams.ny,n);

for i = 1 : n
   q  = rand(dmodel.NB,1);
   dq = rand(dmodel.NB,1);
   
   myBNEA = myBNEA.setState(q, dq);
   [y,yc] = myBNEA.simY();
   myBNEA = myBNEA.setY(y);
   
   bnet{1,i}   = myBNEA.bnt.bnet;
   engine{1,i} = myBNEA.bnt.engine;
   sample(:,i) = yc;
end
% bnet = learn_params_em_modified(engine, sample, 10, 1e-4);

% This parameters should be kept small (1e-15)for
% simulated data. For real data is plays the role
% of a regularization factor and should be tuned so
% as to have non-decresing EM steps.

cov_prior_weight = 1e-15;
bnetHat = EM_bnet_learn(bnet, sample, cov_prior_weight);


