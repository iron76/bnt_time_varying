clear all
close all
clc

if ~exist('preprocess.mat', 'file')
   
   data.nsamples  = 1000; %number of samples
   data.plot      = 0;
   data.ini       = 0;    %seconds to be skipped at the start
   data.end       = 195;  %seconds to reach the end of the movement
   data.diff_imu  = 1;    %derivate the angular velocity of the IMUs
   data.diff_q    = 1;    %derivate the angular velocity of the IMUs


   %%strucutre from files
   data.path        = '/Users/iron/Desktop/iron/myTex/2015-01-rss/data';
   data.parts       = {'inertial'                                 , 'left_arm_accelerometers', 'left_foot_inertial'        , 'left_hand_inertial'                  , 'right_arm_accelerometers', 'right_foot_inertial'       , 'right_hand_inertial'                 , 'torso_accelerometers'                    , 'l_arm_ft_sensor:o', 'r_arm_ft_sensor:o', 'l_leg_ft_sensor:o', 'r_leg_ft_sensor:o', 'l_foot_ft_sensor:o'        , 'r_foot_ft_sensor:o'        , 'head'      , 'left_arm'  , 'right_arm' , 'left_leg'  , 'right_leg' , 'torso'     };
   data.labels      = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_imu'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_imu'                              , 'to_acc'                                  , 'la_fts'           , 'ra_fts'           , 'll_fts'           , 'rl_fts'           , 'lf_fts'                    , 'rf_fts'                    , 'h'         , 'la'        , 'ra'        , 'll'        , 'rl'        , 'to'        };
   data.ndof        = {12                                         , 18                       , 3                           , 6                                     , 18                        , 3                           , 6                                     , 12                                        , 6                  , 6                  , 6                  , 6                  , 6                           , 6                           ,  6          , 16          , 16          , 6           , 6           , 3           };
   data.index       = {'4:9'                                      , '1:3'                    , '1:3'                       , '1:6'                                 , '1:3'                     , '1:3'                       , '1:6'                                 , '1:3'                                     , '1:6'              , '1:6'              , '1:6'              , '1:6'              , '1:6'                       , '1:6'                       ,  '1:6'      , '1:16'      , '1:16'      , '1:6'       , '1:6'       , '1:3'       };
   data.type        = {''                                         , 'analog:o'               , 'analog:o'                  , 'analog:o'                            , 'analog:o'                , 'analog:o'                  , 'analog:o'                            , 'analog:o'                                , ''                 , ''                 , ''                 , ''                 , ''                          , ''                          , 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o'};
   data.visualize   = {1*data.plot                                , 1*data.plot               , 1*data.plot                , 1*data.plot                           , 1*data.plot               , 1*data.plot                 , 1*data.plot                           , 1*data.plot                               , 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot                 , 1*data.plot                 , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot };
   data = loadData(data);
   %%strucutre for urdf
   sens.parts       = {'chest+torso+neck_1+neck_2+head+imu_frame' , 'l_upper_arm+l_arm'      , 'l_upper_foot+l_foot+l_sole', 'l_forearm+l_wrist_1+l_hand+l_gripper', 'r_upper_arm+r_arm'       , 'r_upper_foot+r_foot+r_sole', 'r_forearm+r_wrist_1+r_hand+r_gripper', 'chest+torso+neck_1+neck_2+head+imu_frame', 'l_upper_arm+l_arm', 'r_upper_arm+r_arm', 'l_thigh'          , 'r_thigh'          , 'l_upper_foot+l_foot+l_sole', 'r_upper_foot+r_foot+r_sole'};
   sens.labels      = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_imu'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_imu'                              , 'to_acc'                                  , 'la_fts'           , 'ra_fts'           , 'll_fts'           , 'rl_fts'           , 'lf_fts'                    , 'rf_fts'                    };
   sens.ndof        = {6                                          , 3                        , 3                           , 6                                     , 3                         , 3                           , 6                                     , 3                                         , 6                  , 6                  , 6                  , 6                  , 6                           , 6                           };
   sens.type        = {''                                         , 'analog:o'               , 'analog:o'                  , 'analog:o'                            , 'analog:o'                , 'analog:o'                  , 'analog:o'                            , 'analog:o'                                , ''                 , ''                 , ''                 , ''                 , ''                          , ''                          };

   %%
   % for i = 1 : length(data.time)
   %    iCubVisualize(data.q(:,i), R)
   % end
   
   run('iCub.m')
   dmodel  = iCub_dmodel;
   ymodel  = iCubSens(dmodel, sens);
   
   dmodel  = autoTreeStochastic(dmodel);
   ymodel  = iCubSensStochastic(ymodel);
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
      if mod(i-1,100) == 0
         fprintf('Processing %d %% of the dataset\n', round(i/length(data.time)*100));
      end
   end
   save preprocess.mat
end
%% plot results
load preprocess.mat
close all
py = [0; cumsum(cell2mat(myPNEA.IDsens.sensorsParams.sizes))];
% for i = 1 : myPNEA.IDsens.sensorsParams.ny - dmodel.NB
%    figure
%    J = myPNEA.IDsens.sensorsParams.sizes{i};
%    for j = 1 : J/3
%       subplot([num2str(J/3) '1' num2str(j)])
%       I = py(i)+1+(j-1)*3 : py(i)+3*j;
%       plot(data.time, y(I,:))
%       title(strrep(myPNEA.IDsens.sensorsParams.labels{i}, '_', '~'));
%    end
% end

label_to_plot = {'imu', 'la_acc', 'lf_acc', 'lh_imu', 'ra_acc', 'rf_acc', 'rh_imu', 'to_acc', 'la_fts', 'ra_fts', 'll_fts', 'rl_fts', 'lf_fts', 'rf_fts'};
for l = 1 : length(label_to_plot)
   for i = 1 : length(data.parts)
      if strcmp(data.labels{i}, label_to_plot{l})
         t    = ['time_' data.labels{i}];
         ys   = ['ys_' data.labels{i}];
         
         figure
         J = length(eval(data.index{i}));
         for j = 1 : J/3
            subplot([num2str(J/3) '1' num2str(j)])
            I = 1+(j-1)*3 : 3*j;
            eval(['plot(data.time,data.' ys '(I,:), ''--'' )' ]);
            title(strrep(['y_{' data.labels{i} '}'], '_', '~'))
         end
      end
   end
   
   for k = 1 : myPNEA.IDsens.sensorsParams.ny - dmodel.NB
      if strcmp(myPNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
         figure
         J = myPNEA.IDsens.sensorsParams.sizes{k};
         for j = 1 : J/3
            subplot([num2str(J/3) '1' num2str(j)])
            I = py(k)+1+(j-1)*3 : py(k)+3*j;
            plot(data.time, y(I,:))
            title(strrep(myPNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
         end
      end
   end
end

