clear all
close all
clc

% if ~exist('preprocess.mat', 'file')
   
   data.nsamples  = 100; %number of samples
   data.plot      = 0;
   data.ini       = 50;   %seconds to be skipped at the start
   data.end       = 100;  %seconds to reach the end of the movement
   data.diff_imu  = 1;    %derivate the angular velocity of the IMUs
   data.diff_q    = 1;    %derivate the angular velocity of the IMUs

subjectID = 1;
trialID = 1;
%    %%strucutre from files
%    data.path        = '/Users/iron/Desktop/iron/myTex/2015-01-rss/data';
%    % data.path        = '/home/pegua/Documents/Papers/2015-01-rss/data';
%    data.parts       = {'inertial'                                 , 'left_arm_accelerometers', 'left_foot_inertial'        , 'left_hand_inertial'                  , 'right_arm_accelerometers', 'right_foot_inertial'       , 'right_hand_inertial'                 , 'torso_accelerometers'                    , 'l_arm_ft_sensor:o', 'r_arm_ft_sensor:o', 'l_leg_ft_sensor:o', 'r_leg_ft_sensor:o', 'l_foot_ft_sensor:o'        , 'r_foot_ft_sensor:o'        , 'head'      , 'left_arm'  , 'right_arm' , 'left_leg'  , 'right_leg' , 'torso'     };
%    data.labels      = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_imu'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_imu'                              , 'to_acc'                                  , 'la_fts'           , 'ra_fts'           , 'll_fts'           , 'rl_fts'           , 'lf_fts'                    , 'rf_fts'                    , 'h'         , 'la'        , 'ra'        , 'll'        , 'rl'        , 'to'        };
%    data.ndof        = {12                                         , 18                       , 3                           , 6                                     , 18                        , 3                           , 6                                     , 12                                        , 6                  , 6                  , 6                  , 6                  , 6                           , 6                           ,  6          , 16          , 16          , 6           , 6           , 3           };
%    data.index       = {'4:9'                                      , '1:3'                    , '1:3'                       , '1:6'                                 , '1:3'                     , '1:3'                       , '1:6'                                 , '1:3'                                     , '1:6'              , '1:6'              , '1:6'              , '1:6'              , '1:6'                       , '1:6'                       ,  '1:6'      , '1:16'      , '1:16'      , '1:6'       , '1:6'       , '1:3'       };
%    data.type        = {''                                         , 'analog:o'               , 'analog:o'                  , 'analog:o'                            , 'analog:o'                , 'analog:o'                  , 'analog:o'                            , 'analog:o'                                , ''                 , ''                 , ''                 , ''                 , ''                          , ''                          , 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o', 'stateExt:o'};
%    data.visualize   = {1*data.plot                                , 1*data.plot               , 1*data.plot                , 1*data.plot                           , 1*data.plot               , 1*data.plot                 , 1*data.plot                           , 1*data.plot                               , 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot        , 1*data.plot                 , 1*data.plot                 , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot , 1*data.plot };
%    data = loadData(data);

    %%strucutre for urdf
    %sens.parts       = {'chest+torso+neck_1+neck_2+head+imu_frame' , 'l_upper_arm+l_arm'      , 'l_upper_foot+l_foot+l_sole', 'l_forearm+l_wrist_1+l_hand+l_gripper', 'r_upper_arm+r_arm'       , 'r_upper_foot+r_foot+r_sole', 'r_forearm+r_wrist_1+r_hand+r_gripper', 'chest+torso+neck_1+neck_2+head+imu_frame', 'l_upper_arm+l_arm', 'r_upper_arm+r_arm', 'l_thigh'          , 'r_thigh'          , 'l_upper_foot+l_foot+l_sole', 'r_upper_foot+r_foot+r_sole'};
    %sens.labels      = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_imu'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_imu'                              , 'to_acc'                                  , 'la_fts'           , 'ra_fts'           , 'll_fts'           , 'rl_fts'           , 'lf_fts'                    , 'rf_fts'                    };
    sens.parts       = {'foot','torso'};
    sens.labels      = {'ftx','imu'};
    sens.ndof        = {6,6};
    %sens.type        = {''                                         , 'analog:o'               , 'analog:o'                  , 'analog:o'                            , 'analog:o'                , 'analog:o'                  , 'analog:o'                            , 'analog:o'                                , ''                 , ''                 , ''                 , ''                 , ''                          , ''                          };
 
 %   sens.ndof        = {6                                          , 3                        , 3                           , 6                                     , 3                         , 3                           , 6                                     , 3                                         , 6                  , 6                  , 6                  , 6                  , 6                           , 6                           };
 %   sens.type        = {''                                         , 'analog:o'               , 'analog:o'                  , 'analog:o'                            , 'analog:o'                , 'analog:o'                  , 'analog:o'                            , 'analog:o'                                , ''                 , ''                 , ''                 , ''                 , ''                          , ''                          };
 
    %%structure from sensors 
    %sens.transform   = {'X_chest_imu'                              , 'X_l_upper_arm_la_acc'   , 'X_l_upper_foot_lf_acc'     , 'X_l_forearm_lh_imu'                  , 'X_r_upper_arm_ra_acc'    , 'X_r_upper_foot_rf_acc'     , 'X_r_forearm_rh_imu'                  , 'X_chest_to_acc'                          , 'X_l_upper_arm_la_fts_force','X_r_upper_arm_ra_fts_force','X_l_thigh_ll_fts_force','X_r_thigh_rl_fts_force','X_l_upper_foot_lf_fts_force'     , 'X_r_upper_foot_rf_fts_force'     }; 


%    %%
%    % for i = 1 : length(data.time)
%    %    iCubVisualize(data.q(:,i), R)
%    % end

   %run('humanThreeLink.m')
   humanThreeLink
   dmodel  = humanThreeLink_dmodel;
   ymodel  = humanThreeLinkSens(dmodel, sens);
   
   dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e4);
   ymodel  = humanThreeLinkSensStochastic(ymodel);
   myModel = model(dmodel);
   mySens  = sensors(ymodel);
   myPNEA  = PNEA(myModel, mySens);

   