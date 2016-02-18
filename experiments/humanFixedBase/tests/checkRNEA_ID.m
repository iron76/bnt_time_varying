
% CHECKRNEA_ID.m computes vector d = [d_1,d_2,...,d_NB] comparing two
% equivalent methods: 
% - METHOD 1: using Featherstone ID --> d
% - METHOD 2: using RNEA class --> d_RNEA

clear 
close all
clc

subjectID = 1;
trialID = 1;

%%=====structure from files
data.path        = './experiments/humanFixedBase/intermediateDataFiles/processedSensorData.mat';
[ data ] = organiseBERDYCompatibleSensorData( data, subjectID, trialID );
close all;

data.parts    = {'leg'         ,'torso'};
data.labels   = {'fts'         ,'imu'  };
data.ndof     = {6             ,6      };
data.index    = {'1:6'         ,'1:6'  };

%%=====structure of sensors for URDF
sens.parts    = {'leg'         ,'torso'};  %force of the forceplate is ingoing into the leg
sens.labels   = {'fts'         ,'imu'  };
sens.ndof     = {6             ,6      };

label_to_plot = {'fts'         ,'imu'  };
   
len = length(data.time);

load(sprintf('./experiments/humanFixedBase/data/humanThreeLinkModelFromURDF_subject%d.mat',subjectID));
     
dmodel  = humanThreeLink_dmodel;                     % deterministic model
ymodel  = humanThreeLinkSens(dmodel, sens);  
   
dmodel  = autoTreeStochastic(dmodel, 1e-5, 1e1);     % probabilistic model for D equation (added Sv and Sw)
ymodel  = humanThreeLinkSensStochastic(ymodel);      % probabilistic model for Y(q,dq) d = y (added Sy)
   
myModel = model(dmodel);
mySens  = sensors(ymodel);  

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
 
clear d_temp;
clear tau_i; 
clear a_i;
clear fB_i;
clear f_i;

%% =====METHOD 2: Computing d using Newton-Euler with RNEA method 
tic;

ymodel_RNEA  = autoSensRNEA(dmodel);
mySens_RNEA  = sensors(ymodel_RNEA);
myRNEA       = RNEA(myModel, mySens_RNEA);
      
y_RNEA_f = zeros(6*dmodel.NB, len);
y_RNEA_ddq = zeros(dmodel.NB, len);
fx = cell(dmodel.NB);
  
%Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
for i = 1 : dmodel.NB
   for t = 1 : len
          fx{i,1} = zeros(6,1); 
          y_RNEA_f(((1:6)+(i-1)*6), t) = [fx{i,1}];
          y_RNEA_ddq(i, t) = [data.ddq(i,t)];
   end
   y_RNEA = [y_RNEA_f ; y_RNEA_ddq];
end


d_RNEA = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
   for i = 1 : len
        myRNEA = myRNEA.setState(data.q(:,i), data.dq(:,i));
        myRNEA = myRNEA.setY(y_RNEA(:,i));
        myRNEA = myRNEA.solveID();
       
        d_RNEA(:,i) = myRNEA.d; 
   end
   
t_RNEA = toc;
disp(['[1st] CPU time for d computation with RNEA method is: ' num2str(t_RNEA) '[sec]']);

%% =====check

if (sum(d-d_RNEA) ~= 0);
    disp('Something wrong with d computation. Check methods.');
    res = 1; 
else
    disp('Methods 1 and 2 are equivalent.')
end
