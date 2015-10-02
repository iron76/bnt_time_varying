clear all
close all
NB        = 2;
nsamples  = 1000;

dmodel   = autoTree(NB);
dmodel   = autoTreeStochastic(dmodel);
dmodel.gravity = [0; -9.81; 0];

ymodel  = autoSensRNEA(dmodel);
ymodel  = autoSensStochastic( ymodel );
mySens  = sensors( ymodel );

dmodel  = autoTreeStochastic(dmodel);
myModel = model(dmodel);

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myMAP = MAP(myModel, mySens);
myMAP = myMAP.setState(q,dq);
myMAP = myMAP.solveID('variance');
Sd    = myMAP.Sd;
Sy    = myMAP.IDsens.sensorsParams.Sy.matrix;

[muD, SD] = computePrior(myMAP);
d     = chol(SD)'*randn(26*dmodel.NB, nsamples) + repmat(muD, 1, nsamples);
n     = chol(Sy)'*randn(ymodel.m, nsamples);  % cov(n') = Sy
y     = myMAP.simY(d) + n;
dhat  = zeros(26*dmodel.NB, nsamples);

tic
for i = 1 : nsamples
   myMAP     = myMAP.setY(y(:,i));
   myMAP     = myMAP.solveID('variance');
   dhat(:,i) = myMAP.d;
end
toc

e     = d-dhat;
muHat = mean(e');    % muHat should be zero 
SdHat = cov(e');     % SdHat should be SdHat

% conversions to obtain an lower trinagular [D; Y]
% with these definitions [Y(:,J); D(I, J)] is lower-triangular
I      = [myMAP.ia; myMAP.ifB; myMAP.iF(end:-1:1, 1); myMAP.itau];
J      = [myMAP.jfx; myMAP.jd2q; myMAP.ja; myMAP.jfB; myMAP.jF(end:-1:1, 1); myMAP.jtau];
K      = myMAP.id;
Iinv(I)= 1:length(I);
Jinv(J)= 1:length(J);
Kinv(K)= 1:length(K);

S    = Sd(Kinv, Kinv);
SHat = SdHat(Kinv, Kinv);

figure
imagesc(S(J, J))
colorbar
figure
subplot(211)
imagesc(S([myMAP.jF(end:-1:1, 1)], [myMAP.jF(end:-1:1, 1)]))
title('Force variance')
colorbar
subplot(212)
imagesc(S([myMAP.ja], [myMAP.ja]))
title('Acceleration variance')
colorbar

figure
imagesc(SHat(J, J))
colorbar
figure
subplot(211)
imagesc(SHat([myMAP.jF(end:-1:1, 1)], [myMAP.jF(end:-1:1, 1)]))
title('Force variance')
colorbar
subplot(212)
imagesc(SHat([myMAP.ja], [myMAP.ja]))
title('Acceleration variance')
colorbar

