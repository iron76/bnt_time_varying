clear all
close all

defineConst2;

dtime = 0.40;                %time discretization [s]
[footInertial, tm] = readDataDumper('../../data/dump2/inertial/data.log');
[~, dfootInertial(:,1)] = sgolayDiff(footInertial(:,7), tm);
[~, dfootInertial(:,2)] = sgolayDiff(footInertial(:,8), tm);
[~, dfootInertial(:,3)] = sgolayDiff(footInertial(:,9), tm);
footInertial(:,7:9) = dfootInertial(:,1:3);
time = tm(200):dtime:tm(end-200);    %time samples [s]
footInertial = interp1(tm, footInertial, time, 'spline', 'extrap');
footInertial(:,4:9) = (Adj_IMU*footInertial(:,4:9)')';

[footRightFT, tm] = readDataDumper('../../data/dump2/rightFootFT/data.log');
footRightFT = interp1(tm, footRightFT, time, 'spline', 'extrap');
footRightFT = (Adj_FT*footRightFT')';

[footRightSkin, tm] = readDataDumper('../../data/dump2/rightFootSkin/data.log');
footRightSkin = interp1(tm, footRightSkin, time, 'spline', 'extrap');

[legRightFT, tm] = readDataDumper('../../data/dump2/rightLegFT/data.log');
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

[legRight_q, tm] = readDataDumper('../../data/dump2/rightLegPos/data.log');
[tm, I, J] = unique(tm);
legRight_q = legRight_q(I, :);
legRight_q = interp1(tm, legRight_q, time, 'spline', 'extrap');
legRight_q(:,1) = -legRight_q(:,1)-q0_offset(1);
legRight_q(:,4) = -legRight_q(:,4)-q0_offset(4);

[legRight_dq, tm] = readDataDumper('../../data/dump2/rightLegVel/data.log');
[tm, I, J] = unique(tm);
legRight_dq = legRight_dq(I, :);
legRight_dq = interp1(tm, legRight_dq, time, 'spline', 'extrap');

[legRight_d2q, tm] = readDataDumper('../../data/dump2/rightLegAcc/data.log');
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
[bnet0, index, ns] = init2DOFnet2(n, k, d, l);

%hidden and observed nodes
i_obs    = [...
    index.edo(1) index.edo(2) index.ep(1)   index.ep(2)   index.ec(1)  ...
    index.ec(2)  index.eu(1)  index.eu(2)   index.ef(1)   index.ef(2)  ...
    index.mdo0   index.mp0    index.md2t(1) index.md2t(2) index.mf(1)  ...  
    index.mu(1)  index.mf(3)  index.mu(3)   index.mdo(2)  index.mp(2)  ];

i_cov_learn    = [...
    index.mdo0    index.mdo(1)   index.mdo(2) ...
    index.mp0     index.mp(1)    index.mp(2)  ...
    index.mf(1)   index.mf(2)    index.mf(3)  ...
    index.mu(1)   index.mu(2)    index.mu(3)  ...
    index.md2t(1) index.md2t(2)  ...
    index.mc(1)   index.mc(2)  ];

i_hidden = setdiff(1:k, i_obs);

for j = 1 : nsamples
    %Useful definitions

    f3 = footRightFT(j,1:3);
    u3 = footRightFT(j,4:6);

    f1 = legRightFT(j,1:3);
    u1 = legRightFT(j,4:6);

    do2 =  footInertial(j,7:9);
    p2  = -footInertial(j,4:6);
    
    %Initial conditions (linear and angular accelerations)
    o0x=0; o0y=0; o0z=0;
    p0x=0; p0y=g; p0z=0;
    do0x=0; do0y=0; do0z=0;
    
    o0  = [o0x; o0y; o0z];
    p0  = [p0x; p0y; p0z];
    do0 = [do0x; do0y; do0z];
    
    t   = legRight_q(j,[1 4])*pi/180;
    dt  = legRight_dq(j,[1 4])*pi/180;
    d2t = legRight_d2q(j,[1 4])*pi/180;
    
    bnet{j} = bnet0;
    bnet{j} = build2DOFbnet(bnet{j}, index, t, dt, d2t, o0);
    
    engine{j} = jtree_inf_engine(bnet{j});
    
    % set conditional nodes

    evidence{index.edo(1) ,j}  = zeros(3,1);
    evidence{index.edo(2) ,j}  = zeros(3,1);
    evidence{index.ep(1)  ,j}  = zeros(3,1);
    evidence{index.ep(2)  ,j}  = zeros(3,1);
    evidence{index.ec(1)  ,j}  = zeros(3,1);
    evidence{index.ec(2)  ,j}  = zeros(3,1);
    evidence{index.eu(1)  ,j}  = zeros(3,1);
    evidence{index.eu(2)  ,j}  = zeros(3,1);
    evidence{index.ef(1)  ,j}  = zeros(3,1);
    evidence{index.ef(2)  ,j}  = zeros(3,1);
     
    evidence{index.mdo0   ,j}  = do0;
    evidence{index.mp0    ,j}  = p0;
    evidence{index.md2t(1),j}  = d2t(1);
    evidence{index.md2t(2),j}  = d2t(2);
    evidence{index.mf(1)  ,j}  = f1';
    evidence{index.mu(1)  ,j}  = u1';
    evidence{index.mf(3)  ,j}  = f3';
    evidence{index.mu(3)  ,j}  = u3';
    evidence{index.mdo(2) ,j}  = do2';
    evidence{index.mp(2)  ,j}  = p2';
    
    for h = 1 : length(i_hidden)
        evidence{i_hidden(h), j} = [];
    end
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
    bnet0std{i} = insertStandardization(bnet{i}, M, S, i_cov_learn, covPriorWeight);
    %DO NOT CHANGE THIS WHEN DOING learn_params_em
    
    engineStd{i} = jtree_inf_engine(bnet0std{i});
    [engineStd{i}, ll(i)] = enter_evidence(engineStd{i}, evidence);
    
    engine{i} = jtree_inf_engine(bnet{i});
    [engine{i}, ll(i)] = enter_evidence(engine{i}, evidence);
end

bnetHat = learn_params_em_modified(engineStd, samplesStd, 50, 1e-4);
bnetHat = removeStandardization(bnetHat{nsamples}, M, S, i_cov_learn, covPriorWeight);

%display learning results
hat_d2t1 = struct(bnetHat.CPD{index.md2t(1)});
Sd2t1 = hat_d2t1.cov;

hat_d2t2 = struct(bnetHat.CPD{index.md2t(2)});
Sd2t2 = hat_d2t2.cov;

hat_f1 = struct(bnetHat.CPD{index.mf(1)});
Sf1 = hat_f1.cov;

hat_f2 = struct(bnetHat.CPD{index.mf(2)});
Sf2 = hat_f2.cov;

hat_f3 = struct(bnetHat.CPD{index.mf(3)});
Sf3 = hat_f3.cov;

hat_u1 = struct(bnetHat.CPD{index.mu(1)});
Su1 = hat_u1.cov;

hat_u2 = struct(bnetHat.CPD{index.mu(2)});
Su2 = hat_u2.cov;

hat_u3 = struct(bnetHat.CPD{index.mu(3)});
Su3 = hat_u3.cov;

hat_do0 = struct(bnetHat.CPD{index.mdo0});
Sdo0 = hat_do0.cov;

hat_do1 = struct(bnetHat.CPD{index.mdo(1)});
Sdo1 = hat_do1.cov;

hat_do2 = struct(bnetHat.CPD{index.mdo(2)});
Sdo2 = hat_do2.cov;

hat_p0 = struct(bnetHat.CPD{index.mp0});
Sp0 = hat_p0.cov;

hat_p1 = struct(bnetHat.CPD{index.mp(1)});
Sp1 = hat_p1.cov;

hat_p2 = struct(bnetHat.CPD{index.mp(2)});
Sp2 = hat_p2.cov;

hat_c1 = struct(bnetHat.CPD{index.mc(1)});
Sc1 = hat_c1.cov;

hat_c2 = struct(bnetHat.CPD{index.mc(2)});
Sc2 = hat_c2.cov;


save cov_hat.mat Sd2t1 Sd2t2 Sf1 Sf2 Sf3 Su1 Su2 Su3 Sdo0 Sdo1 Sdo2 Sp0 Sp1 Sp2 Sc1 Sc2


