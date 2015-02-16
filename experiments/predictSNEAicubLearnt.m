clear all
% close all
load preprocess2.mat
load learn_results.mat

[~, n]  = size(y);
NB      = myModel.modelParams.NB;

%% Compute solution
ymodel  = iCubSensLearnt(ymodel, learn);
mySens  = sensors(ymodel);
myPNEA  = PNEA(myModel, mySens);

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
      eval(['res.Sfx_'  link '(:,:,j) = diag(res.Sd(' ind ',' ind ',j));'])
      %d2q
      ind  = '26 + 26*(i-1) : 26*(i-1) + 26';
      eval(['res.d2q_'  link '(:,j)   =  res.d(' ind '        ,j);'])
      eval(['res.Sd2q_' link '(:,:,j) = res.Sd(' ind ',' ind ',j);'])

   end
end

%% Plot overlapped plots
py = [0; cumsum(cell2mat(myPNEA.IDsens.sensorsParams.sizes))];
label_to_plot = {'imu', 'la_acc', 'lf_acc', 'lh_imu', 'ra_acc', 'rf_acc', 'rh_imu', 'to_acc', 'la_fts', 'ra_fts', 'll_fts', 'rl_fts', 'lf_fts', 'rf_fts'};
for l = 1 : length(label_to_plot)
   for k = 1 : myPNEA.IDsens.sensorsParams.ny - dmodel.NB
      if strcmp(myPNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
         figure
         J = myPNEA.IDsens.sensorsParams.sizes{k};
         I = py(k)+1 : py(k)+J;
         colors = ['r', 'g', 'b'];
         for j = 1 : J
            subplot(2, ceil(J/2), j)
            hold on;
            shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);            
            shadedErrorBar(data.time,  res.y(I(j),:), sqrt(reshape(res.Sy(I(j), I(j), :), 1, n)), {colors(mod(j,3)+1) , 'LineWidth', 2}, 0);
            
            title(strrep(myPNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
         end
      end
   end
end

save predictLearnt.mat


