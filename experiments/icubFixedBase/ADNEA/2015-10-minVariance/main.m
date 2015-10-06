close all
clear all
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
Sd1   = myMAP.Sd;

% conversions to obtain an lower trinagular [D; Y]
% with these definitions [Y(:,J); D(I, J)] is lower-triangular
I      = [myMAP.ia; myMAP.ifB; myMAP.iF(end:-1:1, 1); myMAP.itau];
J      = [myMAP.jfx; myMAP.jd2q; myMAP.ja; myMAP.jfB; myMAP.jF(end:-1:1, 1); myMAP.jtau];
K      = myMAP.id;
Iinv(I)= 1:length(I);
Jinv(J)= 1:length(J);
Kinv(K)= 1:length(K);
P      = Kinv(J);
Pinv(P)= 1:length(P);

iLabels = myMAP.iLabels;
jLabels = myMAP.jLabels;
iIndex = myMAP.iIndex;
jIndex = myMAP.jIndex;

figure
imagesc(Sd1(P, P))
colorbar
rangeS = caxis;
jIndex_new  = Pinv(jIndex);
jLabels_new = jLabels(P);
[~, Oj] = sort(jIndex_new);
% iIndex_new  = Iinv(iIndex);
% iLabels_new = iLabels(I);
% [~, Oi] = sort(iIndex_new);
set(gca, 'XTick', jIndex_new(Oj), 'XTickLabel', jLabels_new(jIndex_new(Oj)))
set(gca, 'YTick', jIndex_new(Oj), 'YTickLabel', jLabels_new(jIndex_new(Oj)))
title('Initial covariance matrix')

[~, jMax] = max(diag(abs(Sd1(Kinv, Kinv))));

if ~isempty(find(jMax==myMAP.jF,1))
   num = num2str(find(jMax==myMAP.jF));
   disp(['Worse estimation is on f' num])
   ymodel = addSens(ymodel, ['f' num], eye(6));
   mySens  = sensors( ymodel );
end
if ~isempty(find(jMax==myMAP.jtau,1))
   num = num2str(find(jMax==myMAP.jtau));
   disp(['Worse estimation is on tau' num]);
   ymodel = addSens(ymodel, ['tau' num], 1);
   mySens = sensors( ymodel );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myMAP = MAP(myModel, mySens);
myMAP = myMAP.setState(q,dq);
myMAP = myMAP.solveID('variance');
Sd2   = myMAP.Sd;

% conversions to obtain an lower trinagular [D; Y]
% with these definitions [Y(:,J); D(I, J)] is lower-triangular
I      = [myMAP.ia; myMAP.ifB; myMAP.iF(end:-1:1, 1); myMAP.itau];
J      = [myMAP.jfx; myMAP.jd2q; myMAP.ja; myMAP.jfB; myMAP.jF(end:-1:1, 1); myMAP.jtau];
K      = myMAP.id;
Iinv(I)= 1:length(I);
Jinv(J)= 1:length(J);
Kinv(K)= 1:length(K);
P      = Kinv(J);
Pinv(P)= 1:length(P);

figure
imagesc(Sd2(P, P)) % Sd2(Kinv(J), Kinv(J)) == S(J,J)
colorbar
caxis(rangeS);
jIndex_new  = Pinv(jIndex);
jLabels_new = jLabels(P);
[~, Oj] = sort(jIndex_new);
set(gca, 'XTick', jIndex_new(Oj), 'XTickLabel', jLabels_new(jIndex_new(Oj)))
set(gca, 'YTick', jIndex_new(Oj), 'YTickLabel', jLabels_new(jIndex_new(Oj)))
title('Improved covariance matrix')

figure
imagesc(Sd1(P, P) - Sd2(P, P)) % Sd2(Kinv(J), Kinv(J)) == S(J,J)
colorbar
caxis(rangeS);
jIndex_new  = Pinv(jIndex);
jLabels_new = jLabels(P);
[~, Oj] = sort(jIndex_new);
set(gca, 'XTick', jIndex_new(Oj), 'XTickLabel', jLabels_new(jIndex_new(Oj)))
set(gca, 'YTick', jIndex_new(Oj), 'YTickLabel', jLabels_new(jIndex_new(Oj)))
title('Difference between initial and improved covariance matrix')


