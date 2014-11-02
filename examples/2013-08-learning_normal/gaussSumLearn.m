clc
clear all

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
nsamples = 1000;

i_obs    = [x1 x3];
i_hidden = setdiff(1:N, i_obs);

evidence = cell(1,N);
engine   = cell(1,nsamples);
bnet     = cell(1,nsamples);

bnet = buildGaussSum([10 10], [0 0 0], [10 1 1]);
%DO NOT CHANGE THIS WHEN DOING learn_params_em
engine = jtree_inf_engine(bnet);
[engine, ll] = enter_evidence(engine, evidence);

randn('seed',0);
% %sample the original network
for i=1:nsamples
  samples(:,i) = sample_bnet(bnet, 'evidence', evidence);
end
 
for j = 1 : length(i_hidden)
    for h = 1 : nsamples
        samples{i_hidden(j), h} = [];
    end
end

%normalize for numerical issues 
%samplesStd = (samples - muStd)/covStd
ns  = [1 1 1];
for j = 1 : length(i_obs)
    [dataStd, muStd{i_obs(j)}, covStd{i_obs(j)}] = standardize(cell2mat(samples(i_obs(j), :)));
    samplesStd(i_obs(j), :) = num2cell(dataStd, [1, nsamples]);
end
for j = 1 : length(i_hidden)
    muStd{i_hidden(j)} = zeros(ns(i_hidden(j)),1);
    covStd{i_hidden(j)} = eye(ns(i_hidden(j)));
    for h = 1 : nsamples
        samplesStd{i_hidden(j), h} = [];
    end
end
%samplesStd = S * samples + M
for j = 1 : N
    M{j} = -muStd{j}./covStd{j};
    S{j} = 1./covStd{j};
end

bnet0 = buildGaussSum([10 10], [0 0 0], [10 1 .1]);
bnet0std = insertStandardization(bnet0, M, S);
%DO NOT CHANGE THIS WHEN DOING learn_params_em
engine0std = jtree_inf_engine(bnet0std);
[engine0std, ll] = enter_evidence(engine0std, evidence);

hat_bnetStd = learn_params_em(engine0std, samplesStd);
hat_bnetStd = removeStandardization(hat_bnetStd, M, S);

hat_x1 = struct(hat_bnetStd.CPD{1});
hat_x1.cov

hat_x2 = struct(hat_bnetStd.CPD{2});
hat_x2.cov

hat_x3 = struct(hat_bnetStd.CPD{3});
hat_x3.cov

%Non STD version
% engine0 = jtree_inf_engine(bnet0);
% [engine0, ll] = enter_evidence(engine0, evidence);
% 
% hat_bnet = learn_params_em(engine0, samples);
% 
% hat_x1 = struct(hat_bnet.CPD{1});
% hat_x1.cov
% 
% hat_x2 = struct(hat_bnet.CPD{2});
% hat_x2.cov
% 
% hat_x3 = struct(hat_bnet.CPD{3});
% hat_x3.cov
