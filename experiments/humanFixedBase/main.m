clear 
close all
clc

subjectID = 1;
trialID = 1;

%%
%%
   %%=====structure from files
   data.path        = './experiments/humanFixedBase/data/processedSensorData.mat';
   sensorFrameExtraction
   [data] = organiseBERDYCompatibleSensorData( data, subjectID, trialID );
   close all;
    
   data.parts    = {'leg'         ,'torso'};
   data.labels   = {'fts'         ,'imu'  };
   data.ndof     = {6             ,6      };
   data.index    = {'1:6'         ,'1:6'  };

   %%=====structure of sensors for URDF
   sens.parts    = {'leg'         ,'torso'};             %force of the forceplate is ingoing into the leg
   sens.labels   = {'fts'         ,'imu'  };  
   sens.ndof     = {6             ,6      };
   
   label_to_plot = {'fts'         ,'imu'  };
   
%% Build models 

   load(sprintf('./experiments/humanFixedBase/data/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
   
   dmodel  = humanThreeLink_dmodel;                     %deterministic model
   ymodel  = humanThreeLinkSens(dmodel, sens);  
   
   dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e4);     % probabilistic model for D equation (added Sv and Sw)
   ymodel  = humanThreeLinkSensStochastic(ymodel);      % probabilistic model for Y(q,dq) d = y (added Sy)
   
   myModel = model(dmodel);
   mySens  = sensors(ymodel);  
   
   myMAP  = MAP(myModel, mySens);
   
   len = length(data.time);

%% ================================ RNEA ==================================
%% Compute vector d = [d_1,d_2,...,d_NB] in 2 ways:
%  1) METHOD 1: using Featherstone ID --> d
%  2) METHOD 2: using RNEA class --> d_RNEA
%     d must be == d_RNEA!
%
%  Note: in test checkRNEA.m there is a comparison of ID with iDynTree.

%% ======METHOD 1: Computing d using Newton-Euler with Featherstone ID

tic;

tau = zeros(size(data.q))';
a = cell (size(data.q))';
fB = cell (size(data.q))';
f = cell (size(data.q))';
fx = zeros (6,1);

d_temp = zeros(26*dmodel.NB,1);
d = zeros (26*dmodel.NB, len);

fext    = cell(1,2);
for i = 1 : dmodel.NB
   fext{i}    = fx;
end

for i = 1:len
    
     [tau_i, a_i, fB_i, f_i] = ID(dmodel, data.q(:,i), data.dq(:,i), data.ddq(:,i), fext);
      tau(i,:) = tau_i;
      a(i,:) = a_i;
      fB(i,:) = fB_i;
      f(i,:) = f_i;  
      
      for j = 1 : dmodel.NB
            d_temp((1:26)+(j-1)*26) = [a_i{j}; fB_i{j}; f_i{j}; tau(i,j); fx; data.ddq(j,i)];
      end
      
      d(:,i) = d_temp;
end

t_ID = toc;
disp(['CPU time for d computation with ID method is: ' num2str(t_ID) '[sec]']);
disp('/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ ')


%% =====METHOD 2: Computing d using Newton-Euler with RNEA method
%  
%    tic;
% 
%    ymodel_RNEA  = autoSensRNEA(dmodel);
%    mySens_RNEA  = sensors(ymodel_RNEA);
%    myRNEA       = RNEA(myModel, mySens_RNEA);
%       
%    y_RNEA_f = zeros(6*dmodel.NB, len);
%    y_RNEA_ddq = zeros(dmodel.NB, len);
%    fx = cell(dmodel.NB);
%    
%    %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
%    for i = 1 : dmodel.NB
%       for t = 1 : len
%          fx{i,1} = zeros(6,1); 
%          y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
%          y_RNEA_ddq(i, t) = [data.ddq(i,t)];
%       end
%       y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
%    end
%   
%    d_RNEA = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
%    for i = 1 : len
%        myRNEA = myRNEA.setState(data.q(:,i), data.dq(:,i));
%        myRNEA = myRNEA.setY(y_RNEA(:,i));
%        myRNEA = myRNEA.solveID();
%       
%        d_RNEA(:,i) = myRNEA.d; 
%    end
% 
%   t_RNEA = toc;
%   disp(['[1st] CPU time for tau computation with RNEA method is: ' num2str(t_RNEA) '[sec]']);


%% =====d check

% if (sum(d-d_RNEA) ~= 0);
%    disp('Something wrong with d computation. Check methods.');
%    res = 1;
% else
%    disp('Methods 1 and 2 are equivalent.')
% end
%%
%% Build data.y anda data.Sy 
 
%=====data.y
data.y  = [];
for i = 1 : length(sens.labels)
   eval(['data.y  = [data.y ; data.ys_' sens.labels{i} '];']);
end

% Add the null external forces ftx = 0
data.y  = [data.y; zeros(6*dmodel.NB, len)];
% Add the d2q measurements
data.y  = [data.y; data.ddq];


%=====data.Sy
data.Sy = [];
for i = 1 : length(myMAP.IDsens.sensorsParams.labels)
   data.Sy = [data.Sy; diag(myMAP.IDsens.sensorsParams.Sy{i})];
end

data.Sy = repmat(data.Sy, 1, data.nsamples);
data.Sy = [data.Sy data.Sy(:,end)];


% ordering data.y in angular-linear form -->WRITE BETTER
data.y_cla (1:3,:) = data.y (4:6,:);
data.y_cla (4:6,:) = data.y (1:3,:);
data.y_cla (7:9,:) = data.y (10:12,:);
data.y_cla (10:12,:) = data.y (7:9,:);
data.y_cla (13:26,:) = data.y (13:26,:);

data.y_cla_ord (1:3,:) = data.y (4:6,:);
data.y_cla_ord (4:6,:) = data.y (1:3,:);
data.y_cla_ord (7:9,:) = data.y (10:12,:);
data.y_cla_ord (10:12,:) = data.y (7:9,:);
data.y_cla_ord (13:26,:) = data.y (13:26,:);

% rotating data for making compatible with Drake importing -->WRITE BETTER
data.y_rot = zeros (26, len);

%len =2
for k = 1:len
R_x = [1 0 0; 0 0 1; 0 -1 0];
data.y_rotMom_temp = R_x * data.y_cla(1:3,k);
data.y_rotFor_temp = R_x * data.y_cla(4:6,k);
data.y_rotAcc_temp = R_x * data.y_cla (10:12,k);

data.y_rotMom(:,k) = data.y_rotMom_temp;
data.y_rotFor(:,k) = data.y_rotFor_temp;
data.y_rotAcc(:,k) = data.y_rotAcc_temp;


data.y_rot(1:3,k) = data.y_rotMom(:,k);
data.y_rot(4:6,k) = data.y_rotFor(:,k);
data.y_rot(10:12,k) = data.y_rotAcc(:,k);
data.y_rot(25:26,k) = data.y_cla(25:26,k);

end


%test claudia: ordering data.y_rot in the form [y_1, y_2, ... ,
%y_obj.IDsens.m] -->WRITE BETTER

data.y_rotOrd = zeros(size(data.y_rot));

data.y_rotOrd(1:6,:)= data.y_rot(1:6,:);
data.y_rotOrd(7:12,:)= data.y_rot(13:18,:);
data.y_rotOrd(13,:)= data.y_rot(25,:);
data.y_rotOrd(14:19,:)= data.y_rot(7:12,:);
data.y_rotOrd(20:25,:)= data.y_rot(19:24,:);
data.y_rotOrd(26,:)= data.y_rot(26,:);

%% ================================ MAP ===================================
%% test claudia --> manually build Y -->WRITE BETTER

R_S_0 = [0 -1 0; -1 0 0;0 0 -1];
r_cla = [0.10;0.220;6.03]; %ipotizzato in cm
Xstar_S_0 = [R_S_0  -(R_S_0*skew(r_cla)); zeros(3) R_S_0];

%building Y with transforms -->WRITE BETTER
Y_cla = zeros (26,52);
Y_cla(1:6,13:18) = Xstar_S_0;
Y_cla(10:12,30:32) = processedSensorData.R_2_imu;
Y_cla(13:18,20:25) = zeros(6);
Y_cla(19:24,46:51) = zeros(6);
Y_cla(25,26) = eye(1);
Y_cla(26,52) = eye(1);

%building Y with transforms using measurements in the form [y_1, y_2, ... ,
%y_obj.IDsens.m] -->WRITE BETTER
Y_cla_ord = zeros (26,52);
Y_cla_ord(1:6,13:18) = Xstar_S_0;
Y_cla_ord(7:12,20:25) = zeros(6);
Y_cla_ord(13,26) = zeros(1);
Y_cla_ord(17:19,30:32) = processedSensorData.R_2_imu;
Y_cla_ord(20:25,46:51) = eye(6);
Y_cla_ord(26,52) = eye(1);

%% Computing MAP method

for i = 1 : len
    
    
   myMAP = myMAP.setState(data.q(:,i), data.dq(:,i));
   myMAP = myMAP.setY(data.y_rotOrd(:,i));
   myMAP = myMAP.setYmatrix(Y_cla_ord);
   myMAP = myMAP.solveID();
   
   %Y = cell2mat(myMAP.IDsens.sensorsParams.Y);
  
     res.d(:,i)    = myMAP.d;
     res.Sd(:,:,i) = full(myMAP.Sd); %full() passing from sparse to double matrix
     res.y(:,i)    = Y_cla * res.d(:,i); %non lo capisco
%     res.Sy(:,:,i) = Y * res.Sd(:,:,i) * Y';
%     res.Y(:,:,i) = Y;
% %    
   if mod(i-1,100) == 0
      fprintf('Processing %d %% of the dataset\n', round(i/len*100));
   end
end
% 
% 
% %====plot of Ymatrix
% imagesc(Y)
% colorbar
% title('Y matrix','FontSize',15);


 %% Plot overlapped plots

% py = [0; cumsum(cell2mat(myMAP.IDsens.sensorsParams.sizes))];
% for l = 1 : length(label_to_plot)
%    for k = 1 : myMAP.IDsens.sensorsParams.ny
%       if strcmp(myMAP.IDsens.sensorsParams.labels{k}, label_to_plot{l})
%          figure
%          J = myMAP.IDsens.sensorsParams.sizes{k};
%          I = py(k)+1 : py(k)+J;
%          colors = ['r', 'g', 'b'];
%          for j = 1 : J
%             subplot(2, ceil(J/2), j)
%             hold on;
%             shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);
%             plot(data.time, y(I(j),:), colors(mod(j,3)+1) , 'LineWidth', 1);
% 
%             title(strrep(myMAP.IDsens.sensorsParams.labels{k}, '_', '~'));
%          end
%       end
%    end
% end
%%
% Plotting covariance matrices (some examples)

% %===Covariance Sigma(d|y)
% figure
% subplot(2,1,1)
% imagesc(res.Sd(:,:,1))
% colorbar
% title('Covariance Sigma(d|y) matrix','FontSize',15);
% 
% %===Covariance Sigma(y)
% subplot(2,1,2)
% imagesc(res.Sy(:,:,1))
% colorbar
% title('Covariance Sigma(y) matrix','FontSize',15);


 %% Rerrange solution
 
% d  = zeros(26*dmodel.NB, len);
% Sd = zeros(26*dmodel.NB, 26*dmodel.NB, len);

for i = 1 : dmodel.NB
   for j = 1 : len
       
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

% save(sprintf('./experiments/humanFixedBase/data/savedBERDYresult_subj%d_trial%d.mat',subjectID,trialID));%,'res','data','myMAP');

%% Comparing RNEA/MAP torques

%load ('resultsFromCheckRNEA.mat');

fig = figure();
axes1 = axes('Parent',fig,'FontSize',16);
box(axes1,'on');
hold(axes1,'on');
grid on;

plot1 = plot(data.time,tau(1:len,1), 'lineWidth',2.5); hold on;
set(plot1,'color',[1 0 0]);
plot2 = plot(data.time,tau(1:len,2), 'lineWidth',2.5); hold on;
set(plot2,'color',[0 0.498039215803146 0]);

plot3 = plot(data.time,res.tau_ankle, 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot3,'color',[1 0 0]);
plot4 = plot(data.time,res.tau_hip, 'lineWidth',1.5,'LineStyle','--'); hold on;
set(plot4,'color',[0 0.498039215803146 0]);

leg = legend('$\tau_{1,RNEA}$','$\tau_{2,RNEA}$','$\tau_{1,MAP}$','$\tau_{2,MAP}$','Location','southeast');
set(leg,'Interpreter','latex');
set(leg,'FontSize',18);
xlabel('Time [s]','FontSize',20);
ylabel('Torque [Nm]','FontSize',20);
axis tight;
grid on;


% %berdyResultSensorTest

%% Simulate output in MAP/RNEA
% 
 for  ind = 1:26
     
        figure;

        %====== Comparison between MAP prediction and actual data
        
        %== simulate output in MAP
        y_pred_MAP = myMAP.simY(res.d);
        
        subplot(2,1,1);
        plot1 = plot(data.time,y_pred_MAP(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
        set(plot1,'color',[1 0 0]);
        plot2 = plot(data.time,data.y_rotOrd(ind,:), 'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0 1]);

        leg = legend('MAP Pred','Input data','Location','northeast');
        %set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        %ylabel('Torque[Nm]','FontSize',20);
        title(sprintf('Figure %d',ind));
        axis tight;
        grid on;
        
        %====== %====== Comparison between RNEA prediction and actual data
        
        %== simulate output in RNEA
        y_pred_RNEA = myMAP.simY(d);
         
        subplot(2,1,2);
        plot1 = plot(data.time,y_pred_RNEA(ind,:), 'lineWidth',1.0, 'LineStyle','--'); hold on;
        set(plot1,'color',[1 0 0]);
        plot2 = plot(data.time,data.y_rot(ind,:), 'lineWidth',1.0); hold on;
        set(plot2,'color',[0 0 1]);
        
        leg = legend('RNEA Pred', 'Input data','Location','northeast');
        %set(leg,'Interpreter','latex');
        set(leg,'FontSize',15);
        xlabel('Time [s]','FontSize',15);
        %ylabel('Torque[Nm]','FontSize',20);
        title(sprintf('Figure %d',ind));
        axis tight;
        grid on;

  end

