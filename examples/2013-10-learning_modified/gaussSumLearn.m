clc
clear all

N = 3;
x1 = 1;  x2 = 2;  x3 = 3;  
nsamples = 200;

i_obs    = [x1 x3];
i_hidden = setdiff(1:N, i_obs);

evidence = cell(1,N);
engine   = cell(1,nsamples);
bnet     = cell(1,nsamples);

for i=1:nsamples
    bnet{i} = buildGaussSum([i i], [1/i 1/i 1/i], [10 1 1]);
    %DO NOT CHANGE THIS WHEN DOING learn_params_em
    engine{i} = jtree_inf_engine(bnet{i});
    [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
end

% randn('seed',0);
% %sample the original network
for i=1:nsamples
  samples(:,i) = sample_bnet(bnet{i}, 'evidence', evidence);
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

%different initial conditions
for i=1:nsamples
    bnet0{i} = buildGaussSum([i i], [1/i 1/i 1/i], [10 1 .1]);
    bnet0std{i} = insertStandardization(bnet0{i}, M, S);
    %DO NOT CHANGE THIS WHEN DOING learn_params_em
    
    engineStd{i} = jtree_inf_engine(bnet0std{i});
    [engineStd{i}, ll(i)] = enter_evidence(engineStd{i}, evidence);
    
    engine{i} = jtree_inf_engine(bnet0{i});
    [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
end

bnetHat = learn_params_em_modified(engineStd, samplesStd, 1000, 1e-5);
bnetHat = removeStandardization(bnetHat{nsamples}, M, S);

hat_x1 = struct(bnetHat.CPD{1});
hat_x1.cov

hat_x2 = struct(bnetHat.CPD{2});
hat_x2.cov

hat_x3 = struct(bnetHat.CPD{3});
hat_x3.cov

