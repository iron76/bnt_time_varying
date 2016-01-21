clear 
close all
clc

subjectID = 1;
trialID = 1;

   data.path        = './experiments/humanFixedBase/data/processedSensorData.mat';

   sens.parts       = {'leg','torso'}; %force of the forceplate is ingoing into the leg
   sens.labels      = {'fts','imu'};  
   sens.ndof        = {6,6};

   load(sprintf('./experiments/humanFixedBase/data/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
   dmodel  = humanThreeLink_dmodel; %deterministic model
  
   ymodel  = humanThreeLinkSens(dmodel, sens);  % sModel, sUnkown-covarianceOfd 
   
   dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e4); % probabilistic model for D equation (added Sv and Sw)
   ymodel  = humanThreeLinkSensStochastic(ymodel); % proabilistic model for Y(q,dq) d = y (added Sy)
   
   myModel = model(dmodel);
   mySens  = sensors(ymodel);

   myMAP  = MAP(myModel, mySens);   

%% plot results 

close all
 
% if(exist(data.path,'file'))
%     load(data.path);
% else    
    sensorFrameExtraction
    [ data ] = organiseBERDYCompatibleSensorData( data, subjectID, trialID );
% end
close all

data.parts =  {'leg','torso'};
data.labels = {'fts','imu'  };
data.ndof =   {6    ,6      };
data.index = {'1:6','1:6'};

label_to_plot = {'fts','imu'};

%% Process raw sensor data and bring it in the desired reference frames

% acc_gain = 1; %5.9855e-04;
% deg_to_rad =  1; %pi/180.0;
% gyro_gain = 1; %deg_to_rad*7.6274e-03;
% 
% for l = 1 : length(label_to_plot)
%    for i = 1 : length(data.parts)
%       if strcmp(data.labels{i}, label_to_plot{l})
%          t    = ['time_' data.labels{i}];
%          ys   = ['ys_' data.labels{i}];
%          J = length(eval(data.index{i}));
%          
%          if( strcmp(data.labels{i},'imu') )
%             eval(['data.ys_' data.labels{i} '(4:6,:) = ' ...
%                   'deg_to_rad*data.ys_' data.labels{i} '(4:6,:);']);
%          end
%        
%          if( strcmp(data.labels{i}(end-2:end),'fts') )
%              eval(['data.ys_' data.labels{i} ' = ' ...
%                    'data.ys_' data.labels{i} ';']);
%          end
%       end
%    end
% end

%% Build data.y anda data.Sy 
 
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


%% Computing MAP method

%n=round(length(data.time));
n=length(data.time);

for i = 1 : n
   myMAP = myMAP.setState(data.q(:,i), data.dq(:,i));
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

imagesc(Y)
colorbar
title('Y matrix','FontSize',15);

%% Rerrange solution

d  = zeros(26*dmodel.NB, n);
Sd = zeros(26*dmodel.NB, 26*dmodel.NB, n);

for i = 1 : dmodel.NB
   for j = 1 : n
       
      link = strrep(myMAP.IDmodel.modelParams.linkname{i}, '+', '_');
      joint = strrep(myMAP.IDmodel.modelParams.jointname{i}, '+', '_');
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
      eval(['res.tau_'  joint '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Stau_' joint '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %fx
      ind  = '20 + 26*(i-1) : 26*(i-1) + 25';
      eval(['res.fx_'   link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sfx_'  link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
      %d2q
      ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
      eval(['res.d2q_'  joint '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sd2q_' joint '(:,:,j) = res.Sd(' ind ',' ind ',j);'])
   
   end
end

save(sprintf('./experiments/humanFixedBase/data/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID));%,'res','data','myMAP');


%% Comparing RNEA/MAP torques

load ('resultsFromCheckRNEA.mat');

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,tau(:,1)', 'lineWidth',2.5); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(:,2)', 'lineWidth',2.5); hold on;
set(plot2,'color',[0 0.498039215803146 0]);
plot3 = plot(data.time,res.tau_ankle', 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot3,'color',[1 0 0]);
plot4 = plot(data.time,res.tau_hip', 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot4,'color',[0 0.498039215803146 0]);

leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{1,MAP}$','$\tau_{2,MAP}$','Location','southeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque [Nm]','FontSize',20);
axis tight;
grid on;


%berdyResultSensorTest

% %% Comparing MAP y-pred/
% 
%  for  ind = 1:26
% 
%         y_pred = myMAP.simY(res.d);
% 
%         fig = figure();
%         axes1 = axes('Parent',fig,'FontSize',16);
%         box(axes1,'on');
%         hold(axes1,'on');
%         grid on;
% 
%         plot1 = plot(data.time,y_pred(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
%         set(plot1,'color',[1 0 0]);
%         plot2 = plot(data.time,data.y(ind,:), 'lineWidth',1.0); hold on;
%         set(plot2,'color',[0 0 1]);
% 
%         leg = legend('Map Pred', 'Actual data','Location','northeast');
%         %set(leg,'Interpreter','latex');
%         set(leg,'FontSize',18);
%         xlabel('Time [s]','FontSize',20);
%         %ylabel('Torque[Nm]','FontSize',20);
%         title(sprintf('Figure %d',ind));
%         axis tight;
%         grid on;
% 
% %     figure();
% %     y_pred = myMAP.simY(res.d);
% % 
% %     plot(y_pred(ind,:)); 
% %     hold on; 
% %     plot(data.y(ind,:), '--');
% %     legend('Map Pred', 'Actual data');
% %     title(sprintf('Figure %d',ind));
%   end
