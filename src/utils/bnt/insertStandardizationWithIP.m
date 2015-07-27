function bnetStd = insertStandardizationWithIP(bnet, muStd, covStd, i_cov_learn, covPriorWeight)

bnetStd = bnet;
for i = 1 : length(muStd)
   if( ~isa(bnet.CPD{i}, 'gaussian_net_wrench_CPD') )
       muY_X = struct(bnet.CPD{i}).mean;
   end
   SY_X = struct(bnet.CPD{i}).cov;
   
   i_prts = bnet.parents{i};
   SY  = diag(covStd{i});
   sY  = covStd{i};
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
        % normal CPD case
        W = struct(bnet.CPD{i}).weights;
        if nargin == 5
            if (~isempty(find(i==i_cov_learn,1)))
                bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
            else
                bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
            end
        elseif nargin == 3
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1);
        end
      else
        % custom CPD for inertial parameter estimation
        twist = struct(bnet.CPD{i}).twist;
        inertial_params = struct(bnet.CPD{i}).inertial_params;
        if nargin == 5
            if (~isempty(find(i==i_cov_learn,1)))
                bnetStd.CPD{i} = gaussian_net_wrench_CPD(bnetStd, i, 'twist', twist , 'wrench_weights', sY, 'acceleration_weights', sX, 'inertial_params', inertial_params, 'cov',  SY*SY_X*SY', 'cov_prior_weight', covPriorWeight);
            else
                bnetStd.CPD{i} = gaussian_net_wrench_CPD(bnetStd, i, 'twist', twist , 'wrench_weights', sY, 'acceleration_weights', sX, 'inertial_params', inertial_params, 'cov',  SY*SY_X*SY', 'clamp_inertial_params', 1, 'clamp_cov', 1);
            end
        elseif nargin == 3
            bnetStd.CPD{i} = gaussian_net_wrench_CPD(bnetStd, i, 'twist', twist , 'wrench_weights', sY, 'acceleration_weights', sX, 'inertial_params', inertial_params, 'cov',  SY*SY_X*SY');
        end
      end
   else
      if( ~isa(bnet.CPD{i}, 'gaussian_net_wrench_CPD') )
        % normal CPD case
        if nargin == 5
            if (~isempty(find(i==i_cov_learn,1)))
                bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
            else
                bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
            end
        elseif nargin == 3
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1);
        end
      else
        % custom CPD for inertial parameter estimation
        % the custom CPD always need a continuous parent, so it should not
        % be possible to reach this conditions 
        assert(false)
      end
   end
end
