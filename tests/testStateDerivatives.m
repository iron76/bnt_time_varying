clear all
close all
clc

NB        = 20;
dmodel    = autoTree(NB);

q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);

