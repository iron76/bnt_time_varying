function [bnetHat, ll] = EM_bnet_learn_IP(bnet, sample, cov_prior_weight, i_obs)

[~,n] = size(sample);

bnetStd    = cell(1,n);
engineStd  = cell(1,n);
samStd     = cell(size(sample));

if nargin < 4
   i_obs   = find(~cellfun(@isempty,(sample(:,1))));
end  

muStd   = cell(size(sample(:,1)));
covStd  = cell(size(sample(:,1)));
M       = cell(size(sample(:,1)));
S       = cell(size(sample(:,1)));

ns = bnet{1,1}.node_sizes;
for j = 1 : length(ns)
    M{j} = zeros(ns(j),1);
    S{j} = ones(ns(j),1);
end

for j = 1 : length(i_obs)
    [dataStd, muStd{i_obs(j)}, covStd{i_obs(j)}] = standardize(cell2mat(sample(i_obs(j), :)));
    samStd(i_obs(j), :) = num2cell(dataStd, [1, n]);
    M{i_obs(j)} = -muStd{i_obs(j)}./covStd{i_obs(j)};
    S{i_obs(j)} = 1./covStd{i_obs(j)};
end


for i=1:n
    bnetStd{i} = insertStandardizationWithIP(bnet{i}, M, S, i_obs, cov_prior_weight);
    %DO NOT CHANGE THIS WHEN DOING learn_params_em
    
    engineStd{i} = jtree_inf_engine(bnetStd{i});
    % [~, ll] = enter_evidence(engineStd{i}, samStd(:, i));
end

[bnetHat, ll] = learn_params_em_modified(engineStd, samStd, 10, 1e-4);
bnetHat       = removeStandardizationWithIP(bnetHat{n}, M, S, i_obs, cov_prior_weight);
