
rng(0)
res = 0;

S_dmodel  = 1e-2;
S_ymodel  = 1e-4;

NB_BNEAIP = 4;
dmodel_BNEAIP = autoTree(NB_BNEAIP);
dmodel_BNEAIP = autoTreeStochastic(dmodel_BNEAIP, S_dmodel);
dmodel_BNEAIP.gravity = [0; -9.81; 0];
ymodel_BNEAIP = basicSensingAndAcc(dmodel_BNEAIP);
ymodel_BNEAIP  = autoSensStochastic(ymodel_BNEAIP, S_ymodel);

res = res || testLearnBNEAIP(dmodel_BNEAIP, ymodel_BNEAIP);

if res ~=0 
   disp('[ERROR] One of the tests failed!')
else
   disp('------------------------------------')
   disp('REPORT: All tests ended successfully!')
   disp('------------------------------------')
end
