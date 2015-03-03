clear all
close all
load preprocess2.mat
load d0.mat

[~, n]  = size(y);
NB      = myModel.modelParams.NB;

Sq  = 0.1 * pi/180;
Sdq = 10 * pi/180;

%% Compute solution
myDNEA    = DNEA(myModel, mySens);

res.d  = zeros(26*NB, n);
res.Sd = zeros(26*NB, 26*NB, n);
res.x  = zeros(2 *NB, n);
res.Sx = zeros(2 *NB, 2* NB, n);

for i = 1 : n
   myDNEA = myDNEA.setState(data.q(:,i), data.dq(:,i));
   myDNEA = myDNEA.setY(data.y(:,i));
   myDNEA = myDNEA.setStateVariance(diag([ones(NB,1)*Sq;  ones(NB,1)*Sdq]));
   myDNEA = myDNEA.setD(d0); 
   myDNEA = myDNEA.solveID();
   
   Y = cell2mat(myDNEA.IDsens.sensorsParams.Y);
   
   res.d(:,i)    = myDNEA.d;
   res.x(:,i)    = myDNEA.x;
   res.Sd(:,:,i) = myDNEA.Sd;
   res.Sx(:,:,i) = myDNEA.Sx;
   
   res.y(:,i)    = Y * [res.d(:,i); res.x(:,i)];
   res.Sy(:,:,i) = Y * blkdiag(res.Sd(:,:,i), res.Sx(:,:,i))* Y';
   if mod(i-1,100) == 0
      fprintf('Processing %d %% of the dataset\n', round(i/n*100));
   end
   d0 = res.d(:,i);
end

%% Rerrange solution
for i = 1 : NB
   for j = 1 : n
      link = strrep(myDNEA.IDmodel.modelParams.linkname{i}, '+', '_');
      di   = ['res.d_'   link  '(:,j)'];
      ind  = '1 + 26*(i-1) : 26*(i-1) + 26';
      eval([di '   = res.d(' ind ',j);'])
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

res.q  = res.x(    1:NB, :);
res.dq = res.x(NB+1:end, :);

%% Plot overlapped plots
py = [0; cumsum(cell2mat(myDNEA.IDsens.sensorsParams.sizes))];
label_to_plot = {'imu'                                      , 'la_acc'                 , 'lf_acc'                    , 'lh_acc'                              , 'ra_acc'                  , 'rf_acc'                    , 'rh_acc'                              , 'to_acc'                                  , 'la_fts'                    , 'ra_fts'                   , 'll_fts'               , 'rl_fts'               , 'lf_fts'                    , 'rf_fts'                     , 'lh_gyr'                               , 'rh_gyr'                              };
for l = 1 : length(label_to_plot)
   for k = 1 : myDNEA.IDsens.sensorsParams.ny - dmodel.NB
      if strcmp(myDNEA.IDsens.sensorsParams.labels{k}, label_to_plot{l})
         figure
         J = myDNEA.IDsens.sensorsParams.sizes{k};
         I = py(k)+1 : py(k)+J;
         colors = ['r', 'g', 'b'];
         for j = 1 : J
            subplot(2, ceil(J/2), j)
            hold on;
            shadedErrorBar(data.time, data.y(I(j),:), sqrt(data.Sy(I(j), :)), {[colors(mod(j,3)+1) '--'] , 'LineWidth', 1}, 0);            
            shadedErrorBar(data.time,  res.y(I(j),:), sqrt(reshape(res.Sy(I(j), I(j), :), 1, n)), {colors(mod(j,3)+1) , 'LineWidth', 2}, 0);
            
            title(strrep(myDNEA.IDsens.sensorsParams.labels{k}, '_', '~'));
         end
      end
   end
end
