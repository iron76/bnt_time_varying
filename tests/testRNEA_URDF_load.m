clear all
close all
clc

% edit next like to the install path of drake
path_to_drake_distro = '/home/naveenoid/Workspace/Simulation/Drake/drake-distro/';
addpath(path_to_drake_distro);
addpath(strcat(path_to_drake_distro,'drake'));
addpath(strcat(path_to_drake_distro,'build/matlab'));

addpath_drake;

% currently only works with doublePendulum but can be expanded later
codyco_path = getenv('YARP_DATA_DIRS');
codyco_robot = 'doublePendulum';%getenv('YARP_ROBOT_NAME');
urdfPath = [codyco_path,filesep,'robots',filesep,codyco_robot];
addpath(urdfPath);


NB = 2;

R = RigidBodyManipulator('model.urdf');
dmodel = R.featherstone;
dmodel.jtype = {'R','R'};

ymodel = autoSensRNEA(dmodel);

rmpath_drake;

rmpath(path_to_drake_distro);
rmpath(strcat(path_to_drake_distro,'drake'));
rmpath(strcat(path_to_drake_distro,'build/matlab'));



q         = rand(dmodel.NB,1);
dq        = rand(dmodel.NB,1);
y         = rand(ymodel.m,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myModel = model(dmodel);
mySens  = sensors(ymodel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myRNEA    = RNEA(myModel, mySens);
myRNEA    = myRNEA.setQ(q);
myRNEA    = myRNEA.setDq(dq);
myRNEA    = myRNEA.setY(y);

if (sum(q-myRNEA.IDstate.q))
   error('Something wrong with the setQ method');
end

if (sum(dq-myRNEA.IDstate.dq))
   error('Something wrong with the setDq method');
end

if (sum(y-myRNEA.IDmeas.y))
   error('Something wrong with the setY method');
end

myRNEA = myRNEA.solveID();
myRNEA.d;

for i = 1 : NB
   fx{i}    = y((1:6)+(i-1)*7,1);
   d2q(i,1) = y(7*i);
end

[tau, a, fB, f] = ID( dmodel, q, dq, d2q, fx);

d = zeros(26*NB, 1);
for i = 1 : NB
   d((1:26)+(i-1)*26, 1) = [a{i}; fB{i}; f{i}; tau(i,1); fx{i}; d2q(i,1)];
end

if (sum(d-myRNEA.d))
   error('Something wrong with the solveID method');
end