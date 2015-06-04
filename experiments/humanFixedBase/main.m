clear all
close all
clc

%if ~exist('processedSensorData.mat', 'file')
   
 %  data.nsamples  = 100; %number of samples
 %  data.plot      = 0;
 %  data.ini       = 50;   %seconds to be skipped at the start
 % data.end       = 100;  %seconds to reach the end of the movement
 % data.diff_imu  = 1;    %derivate the angular velocity of the IMUs
 % data.diff_q    = 1;    %derivate the angular velocity of the IMUs

subjectID = 1;
trialID = 1;
%    %%strucutre from files
data.path        = './experiments/humanFixedBase/processedSensorData.mat';
%    % data.path        = '/home/pegua/Documents/Papers/2015-01-rss/data';
data.parts = {'foot','torso'};
data.labels = {'fts','imu'};
data.ndof = {6,6};
data.index = {'1:6','1:6'};

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
    sens.labels      = {'fts','imu'};
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

   

%% plot results
close all
%sens.transform  = {'X_chest_imu'                              , 'X_l_upper_arm_la_acc'   , 'X_l_upper_foot_lf_acc'     , 'X_l_forearm_lh_imu'                  , 'X_r_upper_arm_ra_acc'    , 'X_r_upper_foot_rf_acc'     , 'X_r_forearm_rh_imu'                  , 'X_chest_to_acc'                          , 'X_l_upper_arm_la_fts_force','X_r_upper_arm_ra_fts_force','X_l_thigh_ll_fts_force','X_r_thigh_rl_fts_force','X_l_upper_foot_lf_fts_force'     , 'X_r_upper_foot_rf_fts_force'     }; 
if(exist(data.path,'file'))
    load(data.path);
else    
    sensorFrameExtraction
    [ data ] = organiseBERDYCompatibleSensorData( data );
end

[ data ] = organiseBERDYCompatibleSensorData( data );
label_to_plot = {'fts','imu'};

%% Process raw sensor data and bring it in the desired reference frames
acc_gain = 1.0;%5.9855e-04;
deg_to_rad = pi/180.0;
gyro_gain = 1.0;%deg_to_rad*7.6274e-03;
for l = 1 : length(label_to_plot)
   for i = 1 : length(data.parts)
      if strcmp(data.labels{i}, label_to_plot{l})
         t    = ['time_' data.labels{i}];
         ys   = ['ys_' data.labels{i}];
         J = length(eval(data.index{i}));
%          if( strcmp(data.labels{i},'lh_imu') || ...
%              strcmp(data.labels{i},'rh_imu') )
%             eval(['data.ys_' data.labels{i} '(1:3,:) = ' ...
%                   'acc_gain*data.ys_' data.labels{i} '(1:3,:);']);
%             eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
%                   'gyro_gain*data.ys_' data.labels{i} '(4:6,:);']);
%          end
         if( strcmp(data.labels{i},'imu') )
            eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
                  'deg_to_rad*data.ys_' data.labels{i} '(4:6,:);']);
         end
%           if( strcmp(data.labels{i}(end-2:end),'acc') )
%              eval(['data.ys_' data.labels{i} '(1:3,:) = ' ...
%                    'acc_gain*data.ys_' data.labels{i} '(1:3,:);']);
%              eval(['data.ys_' data.labels{i} ' = ' ...
%                     sens.transform{i} '(1:3,1:3) * ' 'data.ys_' data.labels{i} ';']);
%           end
%          if( strcmp(data.labels{i}(end-2:end),'imu') )
%              eval(['data.ys_' data.labels{i} ' = ' ...
%                    sens.transform{i} ' * ' 'data.ys_' data.labels{i} ';']);
%          % account for the wrong offset present in the input data                
%          elseif( strcmp(data.labels{i}(end-4:end),'f_fts') )
%              eval(['data.ys_' data.labels{i} '(3,:) = ' ...
%                     'data.ys_' data.labels{i} '(3,:) - 3.9;' ]);
%          end
         if( strcmp(data.labels{i}(end-2:end),'fts') )
             eval(['data.ys_' data.labels{i} ' = ' ...
                   'data.ys_' data.labels{i} ';']);
             %eval(['data.ys_' data.labels{i} ' = ' ...
             %      sens.transform{i} ' * ' 'data.ys_' data.labels{i} ';']);
         end
      end
   end
end

%% Build data.y anda data.Sy from adjusted ys_label
data.y  = [];
for i = 1 : length(sens.labels)
   eval(['data.y  = [data.y ; data.ys_' sens.labels{i} '];']);
end

data.Sy = [];
for i = 1 : length(myPNEA.IDsens.sensorsParams.labels)
   data.Sy = [data.Sy; diag(myPNEA.IDsens.sensorsParams.Sy{i})];
end
data.Sy = repmat(data.Sy, 1, data.nsamples);

% Add the null external forces fx = 0
data.y  = [data.y; zeros(6*dmodel.NB, length(data.time))];
% Add the d2q measurements
data.y  = [data.y; data.d2q];

data.Sy = [data.Sy data.Sy(:,end)];
% %% Plot overlapped plots
% py = [0; cumsum(cell2mat(myPNEA.IDsens.sensorsParams.sizes))];
% for l = 1 : length(label_to_plot)
%    for k = 1 : myPNEA.IDsens.sensorsParams.ny
%       if strcmp(myPNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
%          figure
%          J = myPNEA.IDsens.sensorsParams.sizes{k};
%          I = py(k)+1 : py(k)+J;
%          colors = ['r', 'g', 'b'];
%          for j = 1 : J
%             subplot(2, ceil(J/2), j)
%             hold on;
%             shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);
%             plot(data.time, data.y(I(j),:), colors(mod(j,3)+1) , 'LineWidth', 1);
% 
%             title(strrep(myPNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
%          end
%       end
%    end
% end
% 
% %% Plot separated graphs
% % for l = 1 : length(label_to_plot)
% %    for i = 1 : length(data.parts)
% %       if strcmp(data.labels{i}, label_to_plot{l})
% %          t    = ['time_' data.labels{i}];
% %          ys   = ['ys_' data.labels{i}];
% %          
% %          figure
% %          J = length(eval(data.index{i}));
% %          for j = 1 : J/3
% %             subplot([num2str(J/3) '1' num2str(j)])
% %             I = 1+(j-1)*3 : 3*j;
% %             eval(['plot(data.time,data.' ys '(I,:), ''--'' )' ]);
% %             title(strrep(['y_{' data.labels{i} '}'], '_', '~'))
% %          end
% %       end
% %    end
% %    
% %    for k = 1 : myPNEA.IDsens.sensorsParams.ny - dmodel.NB
% %       if strcmp(myPNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
% %          figure
% %          J = myPNEA.IDsens.sensorsParams.sizes{k};
% %          for j = 1 : J/3
% %             subplot([num2str(J/3) '1' num2str(j)])
% %             I = py(k)+1+(j-1)*3 : py(k)+3*j;
% %             plot(data.time, y(I,:))
% %             title(strrep(myPNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
% %          end
% %       end
% %    end
% % end

NB = 2;
n=round(length(data.time));

d  = zeros(26*NB, n);
Sd = zeros(26*NB, 26*NB, n);
for i = 1 : n
   myPNEA = myPNEA.setState(data.q(:,i), data.dq(:,i));
   myPNEA = myPNEA.setY(data.y(:,i));
   myPNEA = myPNEA.solveID();
   
   Y = cell2mat(myPNEA.IDsens.sensorsParams.Y);
   
   res.d(:,i)    = myPNEA.d;
   res.Sd(:,:,i) = myPNEA.Sd;
   
   res.y(:,i)    = Y * res.d(:,i);
   res.Sy(:,:,i) = Y * res.Sd(:,:,i) * Y';
   if mod(i-1,100) == 0
      fprintf('Processing %d %% of the dataset\n', round(i/n*100));
   end
end

%% Rerrange solution
for i = 1 : NB
   for j = 1 : n
      link = strrep(myPNEA.IDmodel.modelParams.linkname{i}, '+', '_');
      di   = ['res.d_'   link  '(:,j)'];
      ind  = '1 + 26*(i-1) : 26*(i-1) + 26';
      eval([di '   = d(' ind ',j);'])
      %a
      ind  = '1 + 26*(i-1) : 26*(i-1) +  6';
      eval(['res.a_'    link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sa_'   link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %fB
      ind  = '7 + 26*(i-1) : 26*(i-1) + 12';
      eval(['res.fB_'   link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.SfB_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %f
      ind  = '13 + 26*(i-1) : 26*(i-1) + 18';
      eval(['res.f_'    link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sf_'   link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %tau
      ind  = '19 + 26*(i-1) : 26*(i-1) + 19';
      eval(['res.tau_'  link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Stau_' link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %fx
      ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
      eval(['res.fx_'   link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sfx_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %d2q
      ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
      eval(['res.d2q_'  link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sd2q_' link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])

   end
end

save('./experiments/humanFixedBase/savedBERDYresult.mat','res','data','myPNEA');

tempT = data.time;%data.time(1:10:end-3);
figure;plot(tempT,res.tau_foot,'b',tempT,res.tau_leg,'r');
xlabel('Time (sec)');
ylabel('Torques (Nm)');
axis tight; hold on;   
legend('\tau_1 (ankle)','\tau_2 (hip)');
title('Torque Estimates');
