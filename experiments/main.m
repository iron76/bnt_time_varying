clear all
close all
clc

data.nsamples  = 1000;  %number of samples
data.plot      = 0;
data.ini       = 42;    %seconds to be skipped at the start
data.end       = 205;   %seconds to reach the end of the movement
%%
data.path        = '/Users/iron/Desktop/iron/myTex/2015-01-rss/data';
data.parts       = {'inertial' , 'head'      , 'left_arm'  , 'right_arm' , 'left_leg'  , 'right_leg' , 'torso'     , 'left_arm_accelerometers', 'left_foot_inertial', 'left_hand_inertial', 'right_arm_accelerometers', 'right_foot_inertial', 'right_hand_inertial', 'torso_accelerometers'};% , 'l_arm_ft_sensor:o', 'r_arm_ft_sensor:o', 'l_leg_ft_sensor:o', 'r_leg_ft_sensor:o', 'l_foot_ft_sensor:o', 'r_foot_ft_sensor:o'};
data.labels      = {'imu'      , 'h'         , 'la'        , 'ra'        , 'll'        , 'rl'        , 'to'        , 'laAcc'                  , 'lfAcc'             , 'lhImu'             , 'raAcc'                   , 'rfAcc'              , 'rhImu'              , 'toAcc'               };%, 'la_ft'            , 'ra_ft'            , 'll_ft'            , 'rl_ft'            , 'lf_ft'             , 'rf_ft'            };
data.ndof        = {12         ,       6     , 16          , 16          , 6           , 6           , 3           , 18                       , 3                   , 6                   , 18                        , 3                    , 6                    , 12                    };%, 6                  , 6                  , 6                  , 6                  , 6                   , 6                  };
data.type        = {''         , 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'analog:o'               , 'analog:o'          , 'analog:o'          , 'analog:o'                , 'analog:o'           , 'analog:o'           , 'analog:o'            };%, ''                 , ''                 , ''                 , ''                 , ''                  , ''                 };
data.visualize   = {1*data.plot, 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot              , 1*data.plot         , 1*data.plot         , 1*data.plot               , 1*data.plot          , 1*data.plot          , 1*data.plot           };%, 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot         , 1*data.plot        };

data = loadData(data);
%%
% for i = 1 : length(data.time)
%    iCubVisualize(data.q(:,i), R)
% end

run('iCub.m')
dmodel  = iCub_dmodel;
ymodel  = iCubSens(dmodel);

dmodel  = autoTreeStochastic(dmodel);
ymodel  = autoSensStochastic(ymodel);
myModel = model(dmodel);
mySens  = sensors(ymodel);
myPNEA  = PNEA(myModel, mySens);

%% RNEA
ymodel_RNEA  = autoSensRNEA(dmodel);
mySens_RNEA  = sensors(ymodel_RNEA);
myRNEA       = RNEA(myModel, mySens_RNEA);

y      = zeros(ymodel.m, length(data.time));
y_RNEA = zeros(ymodel_RNEA.m, length(data.time));

fx = cell(dmodel.NB);
for i = 1 : dmodel.NB
   for t = 1 : length(data.time)
      fx{i,1} = zeros(6,1);
      y_RNEA(1+7*(i-1):7*i, t) = [fx{i,1}; data.d2q(i,t)];
   end
end

%% simulate the output

for i = 1 : length(data.time)
  myRNEA = myRNEA.setState(data.q(:,i), data.dq(:,i));
  myRNEA = myRNEA.setY(y_RNEA(:,i));
  myRNEA = myRNEA.solveID();
  
  myPNEA = myPNEA.setState(data.q(:,i), data.dq(:,i));
  y(:,i) = myPNEA.simY(myRNEA.d);
end

%% plot results
close all
py = [0; cumsum(cell2mat(myPNEA.IDsens.sensorsParams.sizes))];
for i = 1 : myPNEA.IDsens.sensorsParams.ny - dmodel.NB
   figure
   J = myPNEA.IDsens.sensorsParams.sizes{i};
   for j = 1 : J/3
      subplot([num2str(J/3) '1' num2str(j)])
      I = py(i)+1+(j-1)*3 : py(i)+3*j;
      plot(data.time, y(I,:))
      title(strrep(myPNEA.IDsens.sensorsParams.labels{i}, '_', '~'));
   end
end

