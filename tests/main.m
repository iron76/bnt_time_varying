clear all
close all
clc

res       = 0;
NB        = 30;
dmodel    = autoTree(NB);
ymodel    = autoSensRNEA(dmodel);

dmodel    = autoTreeStochastic(dmodel, 1e-2);
ymodel    = autoSensStochastic(ymodel, 1e-5);

res = res || testSubmatrix;

res = res || testSubmatrixSparse;

res = res || testSNEA(dmodel, ymodel);

res = res || testRNEA(dmodel, ymodel);

res = res || testPNEA(dmodel, ymodel);

res = res || testMRNEA(dmodel, ymodel);

res = res || testMNEA(dmodel, ymodel);

res = res || testBNEA(dmodel, ymodel);

res = res || testLearnBNEA(dmodel, ymodel);

res = res || testDerivatives(dmodel, ymodel);

run('iCub.m')
dmodel = iCub_dmodel;
ymodel = iCubSens(dmodel);

dmodel = autoTreeStochastic(dmodel, 1e-3);
ymodel = autoSensStochastic(ymodel, 1e-2);

res = res || testPNEA(dmodel, ymodel);

res = res || testMNEA(dmodel, ymodel);

res = res || testBNEA(dmodel, ymodel);

res = res || testLearnBNEA(dmodel, ymodel);

res = res || testDerivatives(dmodel, ymodel);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended sucessuflly!')
   disp('------------------------------------')
end
