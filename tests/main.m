clear all
close all
clc

res       = 0;
NB        = 30;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;


dmodel_RNEA   = autoTree(NB);
dmodel_RNEA   = autoTreeStochastic(dmodel_RNEA, S_dmodel);
dmodel_RNEA.gravity = [0; -9.81; 0];


ymodel_RNEA   = autoSensRNEA(dmodel_RNEA);
ymodel_RNEA   = autoSensStochastic(ymodel_RNEA, S_ymodel);

dmodel_SNEA   = dmodel_RNEA;
dmodel_SNEA.gravity = [0; -9.81; 0];

ymodel_SNEA   = autoSensSNEA(dmodel_SNEA, {'fx', 'd2q'});
ymodel_SNEA   = autoSensStochastic(ymodel_SNEA, S_ymodel);

dmodel_DNEA = dmodel_SNEA;
ymodel_DNEA = autoSensSNEA(dmodel_DNEA, {'a', 'fx', 'd2q'});
ymodel_DNEA = autoSensDNEA(dmodel_DNEA, ymodel_DNEA, zeros(dmodel_SNEA.NB,1), zeros(dmodel_SNEA.NB, 1));
ymodel_DNEA = autoSensStochastic(ymodel_DNEA, S_ymodel);

dmodel_ANEA = autoTreeStochasticANEA(dmodel_RNEA, S_dmodel);
dmodel_ANEA.gravity = dmodel_SNEA.gravity;

ymodel_ANEA = autoSensANEA(dmodel_ANEA, 0);
ymodel_ANEA = autoSensStochastic(ymodel_ANEA, S_ymodel);

ymodel_DANEA = autoSensANEA(dmodel_ANEA, 1);
ymodel_DANEA  = autoSensDANEA(dmodel_ANEA, ymodel_DANEA, zeros(NB,1), zeros(NB, 1));
ymodel_DANEA  = autoSensStochastic(ymodel_DANEA, S_ymodel);


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

res = res || testRNEA(dmodel_RNEA, ymodel_RNEA);

res = res || testSubmatrix;

res = res || testSubmatrixSparse;

res = res || testSNEA(dmodel_RNEA, ymodel_RNEA);

res = res || testMRNEA(dmodel_RNEA, ymodel_RNEA);

res = res || testPNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testMNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testLearnBNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testBNEA(dmodel_BNEA, ymodel_BNEA);
 
res = res || testDerivativesD(dmodel_RNEA, ymodel_RNEA, dmodel_DNEA, ymodel_DNEA);

res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

res = res || testCalibration(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

res = res || testANEA(dmodel_RNEA, ymodel_RNEA, dmodel_ANEA, ymodel_ANEA);

res = res || testDANEA(dmodel_ANEA, ymodel_ANEA);

res = res || testDANEACalibration(dmodel_ANEA, ymodel_DANEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended successfully!')
   disp('------------------------------------')
end

S_dmodel  = 1e-1;
S_ymodel  = 1e-3;

run('iCub.m')
dmodel_SNEA = iCub_dmodel;
ymodel_SNEA = iCubSens(dmodel_SNEA);

dmodel_SNEA = autoTreeStochastic(dmodel_SNEA, S_dmodel);
ymodel_SNEA = autoSensStochastic(ymodel_SNEA, S_ymodel);

dmodel_DNEA = dmodel_SNEA;
ymodel_DNEA = autoSensDNEA(dmodel_SNEA, ymodel_SNEA, ones(dmodel_SNEA.NB ,1), ones(dmodel_SNEA.NB, 1));
ymodel_DNEA = autoSensStochastic(ymodel_DNEA, S_ymodel);

res = res || testPNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testMNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testBNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testLearnBNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testDerivativesD(dmodel_SNEA, ymodel_SNEA, dmodel_DNEA, ymodel_DNEA);

res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

res = res || testCalibration(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended successfully!')
   disp('------------------------------------')
end
