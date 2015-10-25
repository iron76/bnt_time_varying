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
   
   load(sprintf('./experiments/humanFixedBase/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
   dmodel  = humanThreeLink_dmodel; %deterministic model
  
   ymodel  = humanThreeLinkSens(dmodel, sens); %creating Ys
   
   dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e4);% probabilistic model for D equation (added Sv and Sw)
   ymodel  = humanThreeLinkSensStochastic(ymodel);% proabilistic model for Y(q,dq) d = y (added Sy)
   %ymodel.Sy_inv 
   
   myModel = model(dmodel);
   mySens  = sensors(ymodel);

   %  myPNEA  = PNEA(myModel, mySens);
   myMAP  = MAP(myModel, mySens);
   

%% plot results
close all
%sens.transform  = {'X_chest_imu'                              , 'X_l_upper_arm_la_acc'   , 'X_l_upper_foot_lf_acc'     , 'X_l_forearm_lh_imu'                  , 'X_r_upper_arm_ra_acc'    , 'X_r_upper_foot_rf_acc'     , 'X_r_forearm_rh_imu'                  , 'X_chest_to_acc'                          , 'X_l_upper_arm_la_fts_force','X_r_upper_arm_ra_fts_force','X_l_thigh_ll_fts_force','X_r_thigh_rl_fts_force','X_l_upper_foot_lf_fts_force'     , 'X_r_upper_foot_rf_fts_force'     }; 
if(exist(data.path,'file'))
    load(data.path);
else    
    sensorFrameExtraction
    [ data ] = organiseBERDYCompatibleSensorData( data, subjectID, trialID );
end

data.parts = {'foot','torso'};
data.labels = {'fts','imu'};
data.ndof = {6,6};
data.index = {'1:6','1:6'};

[ data ] = organiseBERDYCompatibleSensorData( data, subjectID,trialID );
label_to_plot = {'fts','imu'};
% humanThreeLink_dmodel.gravity = [+0;9.81;0];

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
         
         if( strcmp(data.labels{i},'imu') )
            eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
                  'deg_to_rad*data.ys_' data.labels{i} '(4:6,:);']);
         end
       
         if( strcmp(data.labels{i}(end-2:end),'fts') )
             eval(['data.ys_' data.labels{i} ' = ' ...
                   'data.ys_' data.labels{i} ';']);
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
for i = 1 : length(myMAP.IDsens.sensorsParams.labels)
   data.Sy = [data.Sy; diag(myMAP.IDsens.sensorsParams.Sy{i})];
end
data.Sy = repmat(data.Sy, 1, data.nsamples);

% Add the null external forces fx = 0
data.y  = [data.y; zeros(6*dmodel.NB, length(data.time))];
% Add the d2q measurements
data.y  = [data.y; data.d2q];

data.Sy = [data.Sy data.Sy(:,end)];

NB = 2;
n=round(length(data.time));

d  = zeros(26*NB, n);
Sd = zeros(26*NB, 26*NB, n);
for i = 1 : n
   myMAP =  myMAP.setState(data.q(:,i), data.dq(:,i));
   myMAP = myMAP.setY(data.y(:,i));
   
   myMAP = myMAP.solveID();
   
   Y = cell2mat(myMAP.IDsens.sensorsParams.Y);
  
   res.d(:,i)    = myMAP.d;
   res.Sd(:,:,i) = myMAP.Sd;
   
   res.y(:,i)    = Y * res.d(:,i);
   res.Sy(:,:,i) = Y * res.Sd(:,:,i) * Y';
   if mod(i-1,100) == 0
      fprintf('Processing %d %% of the dataset\n', round(i/n*100));
   end
end



%% Rerrange solution
for i = 1 : NB
   for j = 1 : n
      link = strrep(myMAP.IDmodel.modelParams.linkname{i}, '+', '_');
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

save(sprintf('./experiments/humanFixedBase/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID),'res','data','myMAP');


%% comparing tau

load ('resultsFromCheckRNEA.mat');

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,tau(:,1), 'lineWidth',2.5); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(:,2), 'lineWidth',2.5); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
plot3 = plot(data.time,res.tau_foot', 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot3,'color',[0 0 0]);
plot4 = plot(data.time,res.tau_leg', 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot4,'color',[0 0 0]);

leg = legend('$\tau_1$','$\tau_2$','$\tau_{MAP}$','Location','northeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque[Nm]','FontSize',20);
axis tight;
grid on;

%vect1 = tau(:,1) - res.tau_foot';
%vect2 = tau(:,2) - res.tau_leg';

%% test Y

figure();
y_pred = myMAP.simY(res.d);

%  plot(y_pred' - data.y')

ind=1;
plot(y_pred(ind,:)); hold on; plot(data.y(ind,:), '--')
