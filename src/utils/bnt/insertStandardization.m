function bnetStd = insertStandardization(bnet, muStd, covStd, i_cov_learn, covPriorWeight)

if(length(covPriorWeight) == 1)
    covPriorWeight = covPriorWeight * ones(length(muStd));
end

bnetStd = bnet;
for i = 1 : length(muStd) %% going down network node by node
   muY_X = struct(bnet.CPD{i}).mean; % existing mean of node i
   SY_X = struct(bnet.CPD{i}).cov; % existing covariance of node i
   
   i_prts = bnet.parents{i};
   SY  = diag(covStd{i});
   muY = muStd{i};
   if(~isempty(i_prts)) % nodes with no parents
      muX = [];
      sX = [];

      for j = 1:length(i_prts)
         muX = [muX; cell2mat(muStd(i_prts(j)))];
         sX  = [sX; cell2mat(covStd(i_prts(j)))];
      end
      SX = diag(sX);
      
      W = struct(bnet.CPD{i}).weights;
      
      if nargin == 5
         if (~isempty(find(i==i_cov_learn,1))) % selected node for updating (i.e. covariance learning)
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX,...
                'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1,...
                'cov_prior_weight', covPriorWeight{i}); % covariance is not clamped but weights and mean are
         else % nodes that need not be updated
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX,...
                'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1,...
                'clamp_cov', 1); % mean, covariance and weights clamped
         end
      elseif nargin == 3 % learn all nodes
         bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX,...
             'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1);
      end
   else % parents of node are empty 
      if nargin == 5
         if (~isempty(find(i==i_cov_learn,1))) % selected node for updating (i.e. covariance learning)
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY, 'cov',  SY*SY_X*SY',...
                'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight{i});
         else % nodes that need not be updated 
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY, 'cov',  SY*SY_X*SY',...
                'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
         end
      elseif nargin == 3 % learn all nodes
         bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY, 'cov',  SY*SY_X*SY',...
             'clamp_mean', 1, 'clamp_weights', 1);
      end
   end
end
