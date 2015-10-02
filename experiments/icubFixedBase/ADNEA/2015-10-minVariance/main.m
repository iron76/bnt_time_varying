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
Sd1   = myMAP.Sd;

% conversions to obtain an lower trinagular [D; Y]
% with these definitions [Y(:,J); D(I, J)] is lower-triangular
I      = [myMAP.ia; myMAP.ifB; myMAP.iF(end:-1:1, 1); myMAP.itau];
J      = [myMAP.jfx; myMAP.jd2q; myMAP.ja; myMAP.jfB; myMAP.jF(end:-1:1, 1); myMAP.jtau];
K      = myMAP.id;
Iinv(I)= 1:length(I);
Jinv(J)= 1:length(J);
Kinv(K)= 1:length(K);

iIndex  = [1 cumsum(myMAP.iSizes)+1];
iLabels = myMAP.iLabels;
jIndex  = [1 cumsum(myMAP.jSizes)+1];
jLabels = myMAP.jLabels;


S    = Sd1(Kinv, Kinv);
figure
imagesc(Sd1(Kinv(J), Kinv(J))) % Sd1(Kinv(J), Kinv(J)) == S(J,J)
colorbar
rangeS = caxis;
set(gca, 'XTick', [], 'XTickLabel', []);
set(gca, 'YTick', [], 'YTickLabel', []);
setLabels('fx',  0*NB+1 : 6 :  6*NB, 'X')
setLabels('fx',  0*NB+1 : 6 :  6*NB, 'Y')
setLabels('q',   6*NB+1            , 'X')
setLabels('q',   6*NB+1            , 'Y')
setLabels('a',   7*NB+1 : 6 : 13*NB, 'X')
setLabels('a',   7*NB+1 : 6 : 13*NB, 'Y')
setLabels('fB', 13*NB+1 : 6 : 19*NB, 'X')
setLabels('fB', 13*NB+1 : 6 : 19*NB, 'Y')
setLabels('f',  19*NB+1 : 6 : 25*NB, 'X', 1)
setLabels('f',  19*NB+1 : 6 : 25*NB, 'Y', 1)
setLabels('t',  25*NB+1            , 'X')
setLabels('t',  25*NB+1            , 'Y')

% figure
% imagesc(S(J, J))
% colorbar
% rangeS = caxis;
% set(gca, 'XTick', jIndex, 'XTickLabel', jLabels{J})

% figure
% subplot(211)
% imagesc(S([myMAP.jF(end:-1:1, 1)], [myMAP.jF(end:-1:1, 1)]))
% colorbar
% rangeF = caxis;
% subplot(212)
% imagesc(S([myMAP.ja], [myMAP.ja]))
% colorbar
% rangeA = caxis;

dS   = diag(abs(S));
[~, jMax] = max(dS);

if ~isempty(find(jMax==myMAP.jF))
   num = num2str(find(jMax==myMAP.jF));
   disp(['Worse estimation is on f' num])
   ymodel = addSens(ymodel, ['f' num], eye(6));
   mySens  = sensors( ymodel );
end
if ~isempty(find(jMax==myMAP.jtau))
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

S    = Sd2(Kinv, Kinv);
figure
imagesc(Sd2(Kinv(J), Kinv(J))) % Sd2(Kinv(J), Kinv(J)) == S(J,J)
colorbar
caxis(rangeS);

set(gca, 'XTick', [], 'XTickLabel', []);
set(gca, 'YTick', [], 'YTickLabel', []);
setLabels('fx',  0*NB+1 : 6 :  6*NB, 'X')
setLabels('fx',  0*NB+1 : 6 :  6*NB, 'Y')
setLabels('q',   6*NB+1            , 'X')
setLabels('q',   6*NB+1            , 'Y')
setLabels('a',   7*NB+1 : 6 : 13*NB, 'X')
setLabels('a',   7*NB+1 : 6 : 13*NB, 'Y')
setLabels('fB', 13*NB+1 : 6 : 19*NB, 'X')
setLabels('fB', 13*NB+1 : 6 : 19*NB, 'Y')
setLabels('f',  19*NB+1 : 6 : 25*NB, 'X', 1)
setLabels('f',  19*NB+1 : 6 : 25*NB, 'Y', 1)
setLabels('t',  25*NB+1            , 'X')
setLabels('t',  25*NB+1            , 'Y')

% figure
% subplot(211)
% imagesc(S([myMAP.jF(end:-1:1, 1)], [myMAP.jF(end:-1:1, 1)]))
% colorbar
% caxis(rangeF);
% subplot(212)
% imagesc(S([myMAP.ja], [myMAP.ja]))
% colorbar
% caxis(rangeA);






