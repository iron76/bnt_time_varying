clear all
close all
clc

mtbSensorCodes =  {'11B1','11B2', ...
   '11B3','11B4', ...
   '11B5', ...
   '11B6', '11B7', ...
   '11B8', '11B9', ...
   '11B10','11B11', ...
   '11B13','11B12'};

mtbSensorLink = {'r_upper_leg','r_upper_leg', ...
   'r_upper_leg','r_upper_leg', ...
   'r_upper_leg',               ...
   'r_upper_leg','r_upper_leg', ...
   'r_lower_leg','r_lower_leg', ...
   'r_lower_leg','r_lower_leg', ...
   'r_foot','r_foot'};

% generate indeces for the mtb sensors:
nrOfMTBAccs = length(mtbSensorLink);
mtbIndices = {};
for i = 1:nrOfMTBAccs
   mtbIndices{i} = strcat(num2str(2+6*(i-1)+4),':',num2str(2+6*i));
end

mtbSensorFrames = {};
for i = 1:nrOfMTBAccs
   mtbSensorFrames{i} = strcat(mtbSensorLink{i},'_acc_mtb_',mtbSensorCodes{i});
end

mtbSensorLabel = {};
for i = 1:nrOfMTBAccs
   mtbSensorLabel{i} = strcat(mtbSensorCodes{i},'_acc');
end


% some sensor are inverted in the model with respect to how are mounted on
% the real robot
mtbInvertedFrames   =  {true,true, ...
   true,true, ...
   true, ...
   false,false, ...
   true,true,   ...
   true,true,   ...
   false,false};



if ~exist('preprocess.mat', 'file')
   
   data.nsamples  = 2000; %number of samples
   data.plot      = 0;
   data.ini       = 2;   %seconds to be skipped at the start
   data.end       = 22;  %seconds to reach the end of the movement
   data.diff_imu  = 1;    %derivate the angular velocity of the IMUs
   data.diff_q    = 1;    %derivate the angular velocity of the IMUs
   
   
   %%strucutre from files
   data.path        = '/Users/traversaro/src/data/dumperYoga9Sep2915';
   data.parts       = {};
   data.labels      = {};
   data.ndof        = {};
   data.index       = {};
   data.type        = {};
   data.visualize   = {};
   %%strucutre for urdf
   sens.parts       = {};
   sens.labels      = {};
   sens.ndof        = {};
   sens.type        = {};
   sens.transform  = {};
   
   %% add ft sensors
   data = addSensToData(data, 'r_leg_ft_sensor:o' , 'rl_fts'  , 6, '1:6', ''           , 1*data.plot);
   sens = addSensToSens(sens, 'r_upper_leg'       , 'rl_fts'  , 6,        ''           ,'drake_r_upper_leg_X_urdf_r_hip_3');
   data = addSensToData(data, 'r_foot_ft_sensor:o', 'rf_fts'  , 6, '1:6', ''           , 1*data.plot);
   sens = addSensToSens(sens, 'r_foot'       , 'rf_fts'  , 6,        ''           ,'drake_r_foot_X_urdf_r_foot');
   
   % add mtb sensors
   for i = 1:nrOfMTBAccs
      sensorTransformName = strcat('drake_', mtbSensorLink{i},'_X_urdf_',mtbSensorFrames{i});
      data = addSensToData(data, 'right_leg/inertialMTB'    , mtbSensorLabel{i}  , 3, mtbIndices{i}, ''           , 1*data.plot);
      sens = addSensToSens(sens, mtbSensorLink{i} , mtbSensorLabel{i}  , 3,        ''           ,sensorTransformName);
   end
   
   %% add joint measurements
   data = addSensToData(data, 'right_leg'         , 'rl'      , 6, '1:6', 'stateExt:o' , 1*data.plot);
   
   data = loadData(data);
   
   
   
   
   %%
   % for i = 1 : length(data.time)
   %    iCubVisualize(data.q(:,i), R)
   % end
   
   run('iCub.m');
   
   % compute the necessary transforms from URDF
   computeURDFToDrakeTransforms;
   
   dmodel  = iCub_dmodel;
   ymodel  = iCubSens(dmodel, sens);
   dmodel  = autoTreeStochastic(dmodel, 1e-1, 1e4);
   ymodel  = iCubSensStochastic(ymodel);
   
   myModel = model(dmodel);
   mySens  = sensors(ymodel);
   myMAP  = MAP(myModel, mySens);
   
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
      
      y(:,i) = myMAP.simY(myRNEA.d);
      if mod(i-1,100) == 0
         fprintf('Processing %d %% of the dataset\n', round(i/length(data.time)*100));
      end
   end
   save preprocess.mat
end
%% plot results
load preprocess.mat
close all

%label_to_plot = {'rl_fts','rf_fts'}
label_to_plot = [mtbSensorLabel,{'rl_fts','rf_fts'}];
%label_to_plot = {'11B13_acc'};

%% Process raw sensor data and bring it in the desired reference frames
acc_gain = 5.9855e-04;
%acc_gain = 1.0;
deg_to_rad = pi/180.0;
gyro_gain = deg_to_rad*7.6274e-03;
for l = 1 : length(label_to_plot)
   for i = 1 : length(data.parts)
      if strcmp(data.labels{i}, label_to_plot{l})
         t    = ['time_' data.labels{i}];
         ys   = ['ys_' data.labels{i}];
         J = length(eval(data.index{i}));
         if( strcmp(data.labels{i},'lh_imu') || ...
               strcmp(data.labels{i},'rh_imu') )
            eval(['data.ys_' data.labels{i} '(1:3,:) = ' ...
               'acc_gain*data.ys_' data.labels{i} '(1:3,:);']);
            eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
               'gyro_gain*data.ys_' data.labels{i} '(4:6,:);']);
         end
         if( strcmp(data.labels{i},'imu') )
            eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
               'deg_to_rad*data.ys_' data.labels{i} '(4:6,:);']);
         end
         if( strcmp(data.labels{i}(end-2:end),'acc') )
            eval(['data.ys_' data.labels{i} '(1:3,:) = ' ...
               'acc_gain*data.ys_' data.labels{i} '(1:3,:);']);
            eval(['data.ys_' data.labels{i} ' = ' ...
               sens.transform{i} '(1:3,1:3) * ' 'data.ys_' data.labels{i} ';']);
            strcat('correcting ',(data.labels{i}),' measures')
            
         end
         if( strcmp(data.labels{i}(end-2:end),'imu') )
            eval(['data.ys_' data.labels{i} ' = ' ...
               sens.transform{i} ' * ' 'data.ys_' data.labels{i} ';']);
            % account for the wrong offset present in the input data
         end
         if( strcmp(data.labels{i}(end-2:end),'fts') )
            eval(['data.ys_' data.labels{i} ' = -normalToStart(' ...
               sens.transform{i} ') * ' 'data.ys_' data.labels{i} ';']);
         end
      end
   end
end

fprintf('Processed raw sensors\n')

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

% for i = 1 : length(data.time)
%    myMAP = myMAP.setState(data.q(:,i), data.dq(:,i));
%    myMAP = myMAP.setY(data.y(:,i));
%    myMAP = myMAP.solveID();
%    y(:,i) = myMAP.simY(myRNEA.d);
%    if mod(i-1,100) == 0
%       fprintf('Processing %d %% of the dataset\n', round(i/length(data.time)*100));
%    end
% end


%% Plot overlapped plots
py = [0; cumsum(cell2mat(myMAP.IDsens.sensorsParams.sizes))];
for l = 1 : length(label_to_plot)
   for k = 1 : myMAP.IDsens.sensorsParams.ny
      if strcmp(myMAP.IDsens.sensorsParams.labels{k}, label_to_plot{l})
         figure
         J = myMAP.IDsens.sensorsParams.sizes{k};
         I = py(k)+1 : py(k)+J;
         colors = ['r', 'g', 'b'];
         for j = 1 : J
            subplot(2, ceil(J/2), j)
            hold on;
            shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);
            plot(data.time, y(I(j),:), colors(mod(j,3)+1) , 'LineWidth', 1);
            
            title(strcat(strrep(myMAP.IDsens.sensorsParams.labels{k}, '_', '~'),num2str(j)));
         end
      end
   end
end

%% Plot state
% figure
% plot(data.time,data.q')
% title('Joint positions')
% figure
% plot(data.time,data.dq')
% title('Joint velocities')
% figure
% plot(data.time,data.d2q')
% title('Joint accelerations')
%

%% Plot separated graphs
% for l = 1 : length(label_to_plot)
%    for i = 1 : length(data.parts)
%       if strcmp(data.labels{i}, label_to_plot{l})
%          t    = ['time_' data.labels{i}];
%          ys   = ['ys_' data.labels{i}];
%
%          figure
%          J = length(eval(data.index{i}));
%          for j = 1 : J/3
%             subplot([num2str(J/3) '1' num2str(j)])
%             I = 1+(j-1)*3 : 3*j;
%             eval(['plot(data.time,data.' ys '(I,:), ''--'' )' ]);
%             title(strrep(['y_{' data.labels{i} '}'], '_', '~'))
%          end
%       end
%    end
%
%    for k = 1 : myMAP.IDsens.sensorsParams.ny - dmodel.NB
%       if strcmp(myMAP.IDsens.sensorsParams.labels{k}, label_to_plot{l})
%          figure
%          J = myMAP.IDsens.sensorsParams.sizes{k};
%          for j = 1 : J/3
%             subplot([num2str(J/3) '1' num2str(j)])
%             I = py(k)+1+(j-1)*3 : py(k)+3*j;
%             plot(data.time, y(I,:))
%             title(strrep(myMAP.IDsens.sensorsParams.labels{k}, '_', '~'));
%          end
%       end
%    end
% end

% save preprocess2.mat
