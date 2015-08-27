function bnet = removeStandardizationWithIP(bnetStd, muStd, covStd, i_cov_learn, covPriorWeight)

bnet = bnetStd;
for i = 1 : length(muStd)
   SY_X = struct(bnetStd.CPD{i}).cov;
   
   i_prts = bnetStd.parents{i};
   SY  = diag(covStd{i});
   sY  = diag(covStd{i});
   muY = muStd{i};
   if(~isempty(i_prts))
      muX = [];
      sX = [];
      for j = 1:length(i_prts)
         muX = [muX; cell2mat(muStd(i_prts(j)))];
         sX  = [sX; cell2mat(covStd(i_prts(j)))];
      end
      SX = diag(sX);
      
      if( ~isa(bnet.CPD{i}, 'gaussian_net_wrench_CPD') )
        muY_X = struct(bnetStd.CPD{i}).mean;
        W = struct(bnetStd.CPD{i}).weights;
        if nargin == 5
            if (~isempty(find(i==i_cov_learn,1)))
                bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY) + W*SX^(-1)*muX, 'cov',  SY^(-1)*SY_X*SY^(-1)', 'weights', SY^(-1)*W*SX, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
            else
                bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY) + W*SX^(-1)*muX, 'cov',  SY^(-1)*SY_X*SY^(-1)', 'weights', SY^(-1)*W*SX, 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
            end
        elseif nargin == 3
            bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY) + W*SX^(-1)*muX, 'cov',  SY^(-1)*SY_X*SY^(-1)', 'weights', SY^(-1)*W*SX, 'clamp_mean', 1, 'clamp_weights', 1);
        end
      else
        % custom CPD for inertial parameters estimation
        twist = struct(bnetStd.CPD{i}).twist;
        inertial_params = struct(bnetStd.CPD{i}).inertial_params;
        if nargin == 5
            if (~isempty(find(i==i_cov_learn,1)))
                bnet.CPD{i} = gaussian_net_wrench_CPD(bnet, i, 'twist', twist , 'force_weights', ones(6,1), 'acceleration_weights', ones(6,1), 'inertial_params', inertial_params, 'cov', SY^(-1)*SY_X*SY^(-1)','cov_prior_weight', covPriorWeight);
            else
                bnet.CPD{i} = gaussian_net_wrench_CPD(bnet, i, 'twist', twist , 'force_weights', ones(6,1), 'acceleration_weights', ones(6,1), 'inertial_params', inertial_params, 'cov',  SY^(-1)*SY_X*SY^(-1)',  'inertial_params', inertial_params, 'clamp_cov', 1);
            end
        elseif nargin == 3
            bnet.CPD{i} = gaussian_net_wrench_CPD(bnet, i, 'twist', twist , 'force_weights', ones(6,1), 'acceleration_weights', ones(6,1), 'inertial_params', inertial_params,  'cov',  SY^(-1)*SY_X*SY^(-1)');
        end   
      end
   else
      if( isa(bnet.CPD{i}, 'gaussian_net_wrench_CPD') )
        % custom CPD for inertial parameter estimation
        % the custom CPD always need a continuous parent, so it should not
        % be possible to reach this conditions 
        assert(false)
      end
      muY_X = struct(bnetStd.CPD{i}).mean;
      if nargin == 5
         if (~isempty(find(i==i_cov_learn,1)))
            bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY)                , 'cov',  SY^(-1)*SY_X*SY^(-1)', 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
         else
            bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY)                , 'cov',  SY^(-1)*SY_X*SY^(-1)', 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
         end
      elseif nargin == 3
         bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY)                , 'cov',  SY^(-1)*SY_X*SY^(-1)', 'clamp_mean', 1, 'clamp_weights', 1);
      end
   end
end