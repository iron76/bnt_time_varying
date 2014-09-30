function bnetStd = insertStandardization(bnet, i_cov_learn, muStd, covStd, covPriorWeight)

bnetStd = bnet;
for i = 1 : length(muStd)
    muY_X = struct(bnet.CPD{i}).mean;
    SY_X = struct(bnet.CPD{i}).cov;
    
    i_prts = bnet.parents{i};
    SY  = diag(covStd{i});
    muY = muStd{i};
    if(~isempty(i_prts))
        muX = [];
        sX = [];
        for j = 1:length(i_prts)
            muX = [muX; cell2mat(muStd(i_prts(j)))];
            sX  = [sX; cell2mat(covStd(i_prts(j)))];
        end
        SX = diag(sX);
        
        W = struct(bnet.CPD{i}).weights;
        
        if (~isempty(find(i==i_cov_learn,1)))
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
        else
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY - SY*W*SX^(-1)*muX, 'cov',  SY*SY_X*SY', 'weights', SY*W*SX^(-1), 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
        end
    else
        if (~isempty(find(i==i_cov_learn,1)))
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1, 'cov_prior_weight', covPriorWeight);
        else
            bnetStd.CPD{i} = gaussian_CPD(bnetStd, i, 'mean', SY*muY_X + muY                   , 'cov',  SY*SY_X*SY', 'clamp_mean', 1, 'clamp_weights', 1, 'clamp_cov', 1);
        end
    end
end
