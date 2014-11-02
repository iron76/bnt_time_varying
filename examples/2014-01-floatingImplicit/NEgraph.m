clear all
close all

defineConst;

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
legRight_dq = -legRight_dq(I, :);
legRight_dq = interp1(tm, legRight_dq, time, 'spline', 'extrap');

[legRight_d2q, tm] = readDataDumper('../../data/dump2/rightLegAcc/data.log');
[tm, I, J] = unique(tm);
legRight_d2q = -legRight_d2q(I, :);
legRight_d2q = interp1(tm, legRight_d2q, time, 'spline', 'extrap');

footInertialHat = zeros(size(footInertial));
legRightFTHat   = zeros(size(legRightFT));

%t0 = 0
time = time - time(1);


bnet = cell(length(time), 1);
[bnet0, index] = init4DOFnet(n, k, d);

for j = 1 : length(time)
    %Useful definitions

    f5 = footRightFT(j,1:3);
    u5 = footRightFT(j,4:6);

    f3 = legRightFT(j,1:3);
    u3 = legRightFT(j,4:6);
    
    do4 =  footInertial(j,7:9);
    p4  = -footInertial(j,4:6);
    
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
    [bnet{j}, R] = build4DOFbnet(bnet{j}, index, t, dt, o0);
    
    engine = jtree_inf_engine(bnet{j});
    
    % set conditional nodes
    evidence = cell(1,k);
     
    evidence{index.mdo0}     = do0;
    evidence{index.mp0}      = p0;
    evidence{index.md2t(1)}  = 0;
    evidence{index.md2t(2)}  = 0;
    evidence{index.md2t(3)}  = d2t(3);
    evidence{index.md2t(4)}  = d2t(4);
    evidence{index.mf(2)}    = R(1:3,1:3,3)*f3';
    evidence{index.mu(2)}    = R(1:3,1:3,3)*u3';
    evidence{index.mf(5)}    = f5';
    evidence{index.mu(5)}    = u5';
    evidence{index.mdo(4)}   = do4';
    evidence{index.mp(4)}    = p4';
    evidence{index.fbb}      = [0 0 0]';
    evidence{index.ubb}      = [0 0 0]';
    
    [engine, ll] = enter_evidence(engine, evidence);
    
    
    marg = marginal_nodes(engine, index.p(4));
    footInertialHat(1:3,j) = reshape(marg.mu,3,1);
    footInertialSigma(1:3,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.do(4));
    footInertialHat(4:6,j) = reshape(marg.mu,3,1);
    footInertialSigma(4:6,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.f(2));
    legRightFTHat(1:3,j) = R(1:3,1:3,3)'*reshape(marg.mu,3,1);
    legRightFTSigma(1:3,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.u(2));
    legRightFTHat(4:6,j) = R(1:3,1:3,3)'*reshape(marg.mu,3,1);
    legRightFTSigma(4:6,j) = diag(reshape(marg.Sigma,3,3));

    margd2t1 = marginal_nodes(engine, index.d2t(1));
    margd2t2 = marginal_nodes(engine, index.d2t(2));
    margd2t3 = marginal_nodes(engine, index.d2t(3));
    margd2t4 = marginal_nodes(engine, index.d2t(4));
    legRightd2tHat(1:4,j) = [margd2t1.mu; margd2t2.mu; margd2t3.mu; margd2t4.mu];
    legRightd2tSigma(1:4,j) = [margd2t1.Sigma; margd2t2.Sigma; margd2t3.Sigma; margd2t4.Sigma];
            
    fprintf(1, 'total time is: %% %.1f\r', j/length(time)*100)
end

close all

%inertial sensor prediction
subplot(221)
plotWithSigma(time, footInertialHat(1:3,:)', footInertialSigma(1:3,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plotWithSigma(time, -footInertial(:,4:6), repmat(diag(Sp2)', length(footInertial(:,4)), 1), '')     %d2p_computed = -d2p_meas
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, footInertialHat(4:6,:)', footInertialSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plotWithSigma(time(1:end), footInertial(:,7:9), repmat(diag(Sdo2)', length(footInertial(:,7)), 1), '')
axis(est_ax);

figure
%FT sensor prediction
subplot(221)
plotWithSigma(time, legRightFTHat(1:3,:)', legRightFTSigma(1:3,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plotWithSigma(time, legRightFT(:,1:3), repmat(diag(Sf1)', length(legRightFT(:,1)), 1), '')
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, legRightFTHat(4:6,:)', legRightFTSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plotWithSigma(time, legRightFT(:,4:6), repmat(diag(Su1)', length(legRightFT(:,4)), 1), '')
axis(est_ax);

figure
%d2t sensor prediction
subplot(221)
plotWithSigma(time, legRightd2tHat(1:2,:)', legRightd2tSigma(1:2,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plotWithSigma(time, 0.*legRight_d2q(:,[1 4])*pi/180, repmat([Sd2t1 Sd2t2], length(legRight_d2q(:,[1 4])), 1), '')
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, legRightd2tHat(3:4,:)', legRightd2tSigma(3:4,:)', '')
est_ax = axis;
subplot(224)
plotWithSigma(time, legRight_d2q(:,[1 4])*pi/180, repmat([Sd2t3 Sd2t4], length(legRight_d2q(:,[1 4])), 1), '')
axis(est_ax);
