function bnetStd = insertStandardization(bnet, muStd, covStd)

bnetStd = bnet;
for i = 1 : length(muStd)
    muY_X = struct(bnet.CPD{i}).mean;
    SY_X = struct(bnet.CPD{i}).cov;
    
    i_prts = bnet.parents{i};
    SY  = diag(covStd{i});
    muY = muStd{i};
    if(~isempty(i_prts))
        SX = diag(cell2mat(covStd(i_prts)));
        muX = cell2mat(muStd(i_prts))';
        W = struct(bnet.CPD{i}).weights;
              
        bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 1e-5);
    else
        bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', 1e-5);
    end
end
