clear all
close all
clc

res       = 0;
NB        = 30;
S_dmodel  = 1e-2;
S_ymodel  = 1e-4;


dmodel_RNEA   = autoTree(NB);
dmodel_RNEA   = autoTreeStochastic(dmodel_RNEA);
dmodel_RNEA.gravity = [0; -9.81; 0];


ymodel_RNEA   = autoSensRNEA(dmodel_RNEA);
ymodel_RNEA   = autoSensStochastic(ymodel_RNEA);

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

dmodel_BNEA    = autoTreeStochastic(dmodel_SNEA, 1, 1);
ymodel_BNEA    = autoSensStochastic(ymodel_SNEA, 1);

disp('Running testRNEA')
res = res || testRNEA(dmodel_RNEA, ymodel_RNEA);

disp('Running testSubmatrix')
res = res || testSubmatrix;

disp('Running testSubmatrixSparse')
res = res || testSubmatrixSparse;

disp('Running testRegressors')
res = res || testRegressors;

disp('Running testSNEA')
res = res || testSNEA(dmodel_RNEA);

disp('Running testConditioning')
res = res || testConditioning(dmodel_RNEA, ymodel_RNEA);

disp('Running testMRNEA')
res = res || testMRNEA(dmodel_RNEA, ymodel_RNEA);

disp('Running testPNEA')
res = res || testPNEA(dmodel_SNEA, ymodel_SNEA);

disp('Running testMNEA')
res = res || testMNEA(dmodel_SNEA, ymodel_SNEA);

disp('Running testLearnBNEA')
res = res || testLearnBNEA(dmodel_SNEA, ymodel_SNEA);

disp('Running testBNEA')
res = res || testBNEA(dmodel_BNEA, ymodel_BNEA);
 
disp('Running testDerivativesD')
res = res || testDerivativesD(dmodel_SNEA, ymodel_SNEA);

disp('Running testDerivatives')
res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

disp('Running testDNEA')
res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

disp('Running testCalibration')
res = res || testCalibration(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

disp('Running testANEA')
res = res || testANEA(dmodel_RNEA, ymodel_RNEA, dmodel_ANEA, ymodel_ANEA);

disp('Running testDANEA')
res = res || testDANEA(dmodel_ANEA, ymodel_ANEA);

disp('Running testDANEACalibration')
res = res || testDANEACalibration(dmodel_ANEA, ymodel_DANEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

disp('Running testDerivativesANEA')
res = res || testDerivativesANEA(dmodel_ANEA, ymodel_ANEA);

disp('Running testMAP')
res = res || testMAP(dmodel_SNEA, ymodel_SNEA);

disp('Running testdMAP')
res = res || testdMAP(dmodel_SNEA, ymodel_SNEA);

disp('Running testLU on RNEA')
res = res || testLU_RNEA(dmodel_RNEA);

disp('Running testLU on ABA')
res = res || testLU_ABA(dmodel_RNEA);

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

res = res || testDerivativesD(dmodel_SNEA, ymodel_SNEA);

res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

res = res || testConditioning(dmodel_SNEA, ymodel_SNEA);

res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

res = res || testCalibration(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA, S_dmodel);

res = res || testMAP(dmodel_SNEA, ymodel_SNEA);

res = res || testdMAP(dmodel_SNEA, ymodel_SNEA);

disp('Running testLU on RNEA')
res = res || testLU_RNEA(dmodel_SNEA);

disp('Running testLU on ABA')
res = res || testLU_ABA(dmodel_SNEA);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended successfully!')
   disp('------------------------------------')
end
