clear all
close all
clc

res       = 0;
NB        = 30;
dmodel_RNEA   = autoTree(NB);
dmodel_RNEA   = autoTreeStochastic(dmodel_RNEA, 1e-2);

ymodel_RNEA   = autoSensRNEA(dmodel_RNEA);
ymodel_RNEA   = autoSensStochastic(ymodel_RNEA, 1e-5);

dmodel_SNEA   = dmodel_RNEA;
dmodel_SNEA.gravity = [0; -9.81; 0];

ymodel_SNEA   = autoSensSNEA(dmodel_SNEA);
ymodel_SNEA   = autoSensStochastic(ymodel_SNEA, 1e-5);

dmodel_DNEA = dmodel_SNEA;
ymodel_DNEA = autoSensDNEA(dmodel_SNEA, ymodel_SNEA, ones(dmodel_SNEA.NB ,1), ones(dmodel_SNEA.NB, 1));
ymodel_DNEA = autoSensStochastic(ymodel_DNEA, 1e-5);

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

res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

res = res || testDNEA(dmodel_DNEA, ymodel_DNEA, dmodel_SNEA, ymodel_SNEA);

res = res || testBNEA(dmodel_BNEA, ymodel_BNEA);

if res ~= 0
   return
end

run('iCub.m')
dmodel_SNEA = iCub_dmodel;
ymodel_SNEA = iCubSens(dmodel_SNEA);

dmodel_SNEA = autoTreeStochastic(dmodel_SNEA, 1e-3);
ymodel_SNEA = autoSensStochastic(ymodel_SNEA, 1e-2);

dmodel_DNEA = dmodel_SNEA;
ymodel_DNEA = autoSensDNEA(dmodel_SNEA, ymodel_SNEA, ones(dmodel_SNEA.NB ,1), ones(dmodel_SNEA.NB, 1));
ymodel_DNEA = autoSensStochastic(ymodel_DNEA, 1e-5);

res = res || testPNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testMNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testBNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testLearnBNEA(dmodel_SNEA, ymodel_SNEA);

res = res || testDerivatives(dmodel_DNEA, ymodel_DNEA);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended successfully!')
   disp('------------------------------------')
end
