function bnet = removeStandardization(bnetStd, muStd, covStd)

bnet = bnetStd;
for i = 1 : length(muStd)
    muY_X = struct(bnetStd.CPD{i}).mean;
    SY_X = struct(bnetStd.CPD{i}).cov;
    
    i_prts = bnetStd.parents{i};
    SY  = diag(covStd{i});
    muY = muStd{i};
    if(~isempty(i_prts))
        SX = diag(cell2mat(covStd(i_prts)));
        muX = cell2mat(muStd(i_prts))';
        W = struct(bnetStd.CPD{i}).weights;
              
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY) + W*SX^(-1)*muX, 'cov',  SY^(-1)*SY_X*SY^(-1)', 'weights', SY^(-1)*W*SX, 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 1e-5);
    else
        bnet.CPD{i} = gaussian_CPD(bnet, i, 'mean', SY^(-1)* (muY_X - muY)                , 'cov',  SY^(-1)*SY_X*SY^(-1)', 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 1e-5);
    end
end
