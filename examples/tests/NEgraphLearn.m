run NEgraph

i_obs    = [i_f(1) i_u(1) i_do0 i_p0 i_f(3) i_u(3) i_d2t];
i_hidden = setdiff(1:k, i_obs);

%uncertainties
%sModel  = 1;
%sUknown = 1;
%sKnown  = 1;

bnet0 = mk_bnet(new_dag, new_ns, 'discrete', dnodes, 'observed', i_obs);

bnet0.CPD{i_do0}    = gaussian_CPD(bnet, i_do0,    ...
    'mean', do0,    'cov', 1/100.*sUknown.*eye(d,d)); %, 'clamp_mean', 1, 'clamp_cov', 1);
bnet0.CPD{i_p0}     = gaussian_CPD(bnet, i_p0,     ...
    'mean', p0,     'cov', 10.*sUknown.*eye(d,d));% ,    'clamp_mean', 1, 'clamp_cov', 1);
bnet0.CPD{i_d2t(1)} = gaussian_CPD(bnet, i_d2t(1), ...
    'mean', d2t(1), 'cov', sUknown); %,                  'clamp_mean', 1, 'clamp_cov', 1);
bnet0.CPD{i_d2t(2)} = gaussian_CPD(bnet, i_d2t(2), ...
    'mean', d2t(2), 'cov', sUknown); %,                  'clamp_mean', 1, 'clamp_cov', 1);
bnet0.CPD{i_f(3)}   = gaussian_CPD(bnet, i_f(3),   ...
    'mean', f3,     'cov', sUknown.*eye(d,d)); % ,        'clamp_mean', 1, 'clamp_cov', 1);
bnet0.CPD{i_u(3)}   = gaussian_CPD(bnet, i_u(3),   ...
    'mean', u3,     'cov', 1/10.*sUknown.*eye(d,d)); % , 'clamp_mean', 1, 'clamp_cov', 1);

%Kinematic recursion
for i = 1:2
    if i==1
        Ri=R1;
    else
        Ri=R2;
    end
    if i==1
        o(:,i)  = Ri'*(o0+dt(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        b = Ri'*dt(i)*cross(o0,z0);
        if (i_d2t(i) < i_do0)
            W  = [Ri'*z0 Ri'];
            % Dw = diag([covStd{i_d2t(i)}; covStd{i_do0}]);
        else
            W  = [Ri' Ri'*z0];
            % Dw = diag([covStd{i_do0}; covStd{i_d2t(i)}]);
        end
        
        bnet0.CPD{i_do(i)}  = gaussian_CPD(bnet0, i_do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W,...
            'clamp_weights', 1, 'clamp_mean', 1);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        if (i_p0 < i_do(i))
            W = [Ri' -vec_hat(r0(:,i))];
        else
            W = [-vec_hat(r0(:,i)) Ri'];
        end
        bnet0.CPD{i_p(i)}  = gaussian_CPD(bnet0, i_p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));        
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        if (i_p(i) < i_do(i))
            W = [eye(d) -vec_hat(rc(:,i))];
        else
            W = [-vec_hat(rc(:,i)) eye(d)];
        end        
        bnet0.CPD{i_c(i)}  = gaussian_CPD(bnet0, i_c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
        
    else
        o(:,i)  = Ri'*(o(:,i-1)+dt(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        b = Ri'*dt(i)*cross(o(:,i-1),z0);
        if (i_d2t(i) < i_do0)
            W = [Ri'*z0 Ri'];
        else
            W = [Ri' Ri'*z0];
        end
        bnet0.CPD{i_do(i)}  = gaussian_CPD(bnet0, i_do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        if (i_p(i-1) < i_do(i))
            W = [Ri' -vec_hat(r0(:,i))];
        else
            W = [-vec_hat(r0(:,i)) Ri'];
        end
        bnet0.CPD{i_p(i)}  = gaussian_CPD(bnet0, i_p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));        
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        if (i_p(i) < i_do(i))
            W = [eye(d) -vec_hat(rc(:,i))];
        else
            W = [-vec_hat(rc(:,i)) eye(d)];
        end   
        bnet0.CPD{i_c(i)}  = gaussian_CPD(bnet0, i_c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
    end
        
end

% Dynamic recursion
for i = 2:-1:1
    if i==2
        Ri=R3;
        Ii = I2;
    else
        Ri=R2;
        Ii = I1;
    end
    %% f(:,i)  = Ri'*f(:,i+1) + m(i,1).*c(:,i);
    b = zeros(d,1);
    if (i_f(i+1) < i_c(i))
        W = [Ri eye(d).*m(i,1)];
    else
        W = [eye(d).*m(i,1) Ri];
    end
    bnet0.CPD{i_f(i)}  = gaussian_CPD(bnet0, i_f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
    
    %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
    %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i));
    b = cross(o(:,i), Ii*o(:,i));
    
    Wt(:, :, 1) = Ii;
    Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
    Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
    Wt(:, :, 4) = Ri;
    [y,ind] = sort([i_do(i), i_f(i), i_f(i+1),  i_u(i+1)]);
    W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))]; 
    
    bnet0.CPD{i_u(i)}  = gaussian_CPD(bnet0, i_u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W, ...
            'clamp_weights', 1, 'clamp_mean', 1);
end

engine0 = jtree_inf_engine(bnet0);

evidence = cell(1,k);
[engine0, ll] = enter_evidence(engine0, evidence);

nsamples = 100;
%sample the original network
for i=1:nsamples
  samples(:,i) = sample_bnet(bnet, 'evidence', evidence);
end

for j = 1 : length(i_obs)
    [dataStd, muStd{i_obs(j)}, covStd{i_obs(j)}] = standardize(cell2mat(samples(i_obs(j), :)));
    samplesStd(i_obs(j), :) = num2cell(dataStd, [1, nsamples]);
    
end

for j = 1 : length(i_hidden)
    muStd{i_hidden(j)} = zeros(new_ns(i_hidden(j)),1);
    covStd{i_hidden(j)} = eye(new_ns(i_hidden(j)),1);
    for h = 1 : nsamples
        samplesStd{i_hidden(j), h} = [];
    end
end

bnet2 = learn_params_em(engine0, samplesStd);

D = diag(covStd{i_do0});
do0   = struct(bnet.CPD{i_do0});       mu_do0 = do0.mean;                    cov_do0 = do0.cov;
do0_2 = struct(bnet2.CPD{i_do0});    mu_do0_2 = do0_2.mean + muStd{i_do0}; cov_do0_2 = D*do0_2.cov*D;

do1   = struct(bnet.CPD{i_do(1)});     mu_do1 = do1.mean;     cov_do1 = do1.cov;       W_do1 = do1.weights;
do1_2 = struct(bnet2.CPD{i_do(1)});  mu_do1_2 = do1_2.mean; cov_do1_2 = D*do1_2.cov*D; W_do1_2 = do1_2.weights;

do2   = struct(bnet.CPD{i_do(2)});     mu_do2 = do2.mean;     cov_do2 = do2.cov;
do2_2 = struct(bnet2.CPD{i_do(2)});  mu_do2_2 = do2_2.mean; cov_do2_2 = D*do2_2.cov*D;

D = diag(covStd{i_d2t(1)});
d2t1   = struct(bnet.CPD{i_d2t(1)});     mu_d2t1 = d2t1.mean;     cov_d2t1 = d2t1.cov;
d2t1_2 = struct(bnet2.CPD{i_d2t(1)});  mu_d2t1_2 = d2t1_2.mean; cov_d2t1_2 = D*d2t1_2.cov*D;

D = diag(covStd{i_d2t(2)});
d2t2   = struct(bnet.CPD{i_d2t(2)});     mu_d2t2 = d2t2.mean;     cov_d2t2 = d2t2.cov;
d2t2_2 = struct(bnet2.CPD{i_d2t(2)});  mu_d2t2_2 = d2t2_2.mean; cov_d2t2_2 = D*d2t2_2.cov*D;

D = diag(covStd{i_f(3)});
f3    = struct(bnet.CPD{i_f(3)});       mu_f3 = f3.mean;                    cov_f3 = f3.cov;
f3_2  = struct(bnet2.CPD{i_f(3)});    mu_f3_2 = f3_2.mean + muStd{i_f(3)}; cov_f3_2 = D*f3_2.cov*D;

