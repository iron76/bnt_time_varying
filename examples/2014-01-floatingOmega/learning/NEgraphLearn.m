clear all
close all

defineConst;

dtime = 0.40;                %time discretization [s]
[footInertial, tm] = readDataDumper('../../../data/dump2/inertial/data.log');
[footInertial(:,7), dfootInertial(:,1)] = sgolayDiff(footInertial(:,7), tm);
[footInertial(:,8), dfootInertial(:,2)] = sgolayDiff(footInertial(:,8), tm);
[footInertial(:,9), dfootInertial(:,3)] = sgolayDiff(footInertial(:,9), tm);
footInertial(:,7:9)  = pi/180.*footInertial(:,7:9);
dfootInertial(:,1:3) = pi/180.*(Adj_IMU(1:3, 1:3) * dfootInertial(:,1:3)')';
time = tm(300):dtime:tm(end-300);    %time samples [s]
footInertial  = interp1(tm,  footInertial, time, 'spline', 'extrap');
dfootInertial = interp1(tm, dfootInertial, time, 'spline', 'extrap');
footInertial(:,4:9) = (Adj_IMU*footInertial(:,4:9)')';

[footRightFT, tm] = readDataDumper('../../../data/dump2/rightFootFT/data.log');
footRightFT = interp1(tm, footRightFT, time, 'spline', 'extrap');
footRightFT = (Adj_FT*footRightFT')';

[footRightSkin, tm] = readDataDumper('../../../data/dump2/rightFootSkin/data.log');
footRightSkin = interp1(tm, footRightSkin, time, 'spline', 'extrap');

[legRightFT, tm] = readDataDumper('../../../data/dump2/rightLegFT/data.log');
legRightFT = interp1(tm, legRightFT, time, 'spline', 'extrap');
legRightFT = (Adj_FT*legRightFT')';

if (REMOVE_OFFSETS)
    finit = [0; 0; 0];
    uinit = [0; 0; 0];
    footRightFT = footRightFT + repmat([finit; uinit]'-footRightFT(1,:), length(footRightFT'), 1);
    
    finit = Rz(q0_offset(1)*pi/180)*(m1+m2)*g0;
    uinit = [0; 0; 0];
    legRightFT = legRightFT + repmat([finit; uinit]'-legRightFT(1,:), length(legRightFT'), 1);   
end

[legRight_q, tm] = readDataDumper('../../../data/dump2/rightLegPos/data.log');
[tm, I, J] = unique(tm);
legRight_q = legRight_q(I, :);
legRight_q = interp1(tm, legRight_q, time, 'spline', 'extrap');
legRight_q(:,1) = -legRight_q(:,1)-q0_offset(1);
legRight_q(:,4) = -legRight_q(:,4)-q0_offset(4);

[legRight_dq, tm] = readDataDumper('../../../data/dump2/rightLegVel/data.log');
[tm, I, J] = unique(tm);
legRight_dq = legRight_dq(I, :);
legRight_dq = interp1(tm, legRight_dq, time, 'spline', 'extrap');

[legRight_d2q, tm] = readDataDumper('../../../data/dump2/rightLegAcc/data.log');
[tm, I, J] = unique(tm);
legRight_d2q = legRight_d2q(I, :);
legRight_d2q = interp1(tm, legRight_d2q, time, 'spline', 'extrap');

footInertialHat = zeros(size(footInertial));
legRightFTHat   = zeros(size(legRightFT));

%t0 = 0
time = time - time(1);
nsamples = length(time);

%allocate memory
evidence = cell(k,nsamples);
engine   = cell(1,nsamples);
bnet     = cell(1,nsamples);

%init network (constant part)
[bnet0, index, ns] = init4DOFnet(n, k, d);

%hidden and observed nodes
i_obs    = [index.md2t(3) index.md2t(4) index.mf(2)  index.mu(2)...  
            index.mf(5)  index.mu(5)   index.mo(4)  index.mp(4)];

i_hidden = setdiff(1:k, i_obs);

for j = 1 : nsamples
    %Useful definitions

    f5 = footRightFT(j,1:3)';
    u5 = footRightFT(j,4:6)';

    f3 = legRightFT(j,1:3)';
    u3 = legRightFT(j,4:6)';
    
    o4  =  footInertial(j,7:9)';
    do4 = dfootInertial(j,1:3)';
    p4  = -footInertial(j,4:6)';
    
    %Initial conditions (linear and angular accelerations)
    ox=0;  oy=0;   oz=0;
    px=0;  py=g;  pz=0;
    dox=0; doy=0; doz=0;
    
    ob  = [ox; oy; oz];
    pb  = [px; py; pz];
    dob = [dox; doy; doz];
    
    o0  = Rb0'*ob;
    p0  = Rb0'*pb;
    do0 = Rb0'*dob;
    
    t   = [0 0 legRight_q(j,[1 4])*pi/180];
    dt  = [0 0 legRight_dq(j,[1 4])*pi/180];
    d2t = [0 0 legRight_d2q(j,[1 4])*pi/180];
    
    bnet{j}      = bnet0;
    [bnet{j}, R] = build4DOFbnet(bnet{j}, index, t, dt, o0, do0, p0);
    
    engine{j} = jtree_inf_engine(bnet{j});
     
    evidence{index.fbb    , j} = [0 0 0]';
    evidence{index.ubb    , j} = [0 0 0]';
    evidence{index.fee    , j} = [0 0 0]';
    evidence{index.uee    , j} = [0 0 0]';   
    
    evidence{index.md2t(3), j}  = d2t(3);
    evidence{index.md2t(4), j}  = d2t(4);    
    evidence{index.mf(2)  , j}  = R(1:3,1:3,3)*f3;
    evidence{index.mu(2)  , j}  = R(1:3,1:3,3)*u3;
    evidence{index.mf(5)  , j}  = f5;
    evidence{index.mu(5)  , j}  = u5;
    evidence{index.mo(4) , j}   = o4;
    evidence{index.mp(4)  , j}  = p4;
        
    for h = 1 : length(i_hidden)
        evidence{i_hidden(h), j} = [];
    end
    fprintf(1, 'total time is: %% %.1f\r', j/length(time)*100)
end

for j = 1 : length(i_obs)
    [dataStd, muStd{i_obs(j)}, covStd{i_obs(j)}] = standardize(cell2mat(evidence(i_obs(j), :)));
    samplesStd(i_obs(j), :) = num2cell(dataStd, [1, nsamples]);
end
for j = 1 : length(i_hidden)
    muStd{i_hidden(j)} = zeros(ns(i_hidden(j)),1);
    covStd{i_hidden(j)} = ones(ns(i_hidden(j)), 1);
    for h = 1 : nsamples
        samplesStd{i_hidden(j), h} = [];
    end
end
%samplesStd = S * samples + M
for j = 1 : k
    M{j} = -muStd{j}./covStd{j};
    S{j} = 1./covStd{j};
end

for i=1:nsamples
    bnet0std{i} = insertStandardization(bnet{i}, i_obs, M, S, covPriorWeight);
    %DO NOT CHANGE THIS WHEN DOING learn_params_em
    
    engineStd{i} = jtree_inf_engine(bnet0std{i});
    [engineStd{i}, ll(i)] = enter_evidence(engineStd{i}, evidence);
    
    engine{i} = jtree_inf_engine(bnet{i});
    [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
end

save dataElaboration.mat

bnetHat = learn_params_em_modified(engineStd, samplesStd, 10, 1e-4);
bnetHat = removeStandardization(bnetHat{nsamples}, i_obs, M, S, covPriorWeight);

%display learning results
hat_d2t3 = struct(bnetHat.CPD{index.md2t(3)});
Sd2t3 = hat_d2t3.cov;

hat_d2t4 = struct(bnetHat.CPD{index.md2t(4)});
Sd2t4 = hat_d2t4.cov;

hat_f2 = struct(bnetHat.CPD{index.mf(2)});
Sf2 = hat_f2.cov;

hat_f5 = struct(bnetHat.CPD{index.mf(5)});
Sf5 = hat_f5.cov;

hat_u2 = struct(bnetHat.CPD{index.mu(2)});
Su2 = hat_u2.cov;

hat_u5 = struct(bnetHat.CPD{index.mu(5)});
Su5 = hat_u5.cov;

hat_o4 = struct(bnetHat.CPD{index.mo(4)});
So4 = hat_o4.cov;

hat_p4 = struct(bnetHat.CPD{index.mp(4)});
Sp4 = hat_p4.cov;

save cov_hat.mat ...
    Sf2   Sf5 ...
    Su2   Su5 ...
    So4  Sp4  ...
    Sd2t3 Sd2t4 


