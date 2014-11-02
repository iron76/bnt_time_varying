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


bnet = cell(length(time), 1);
[bnet0, index] = init2DOFnet3(n, k, d);

for j = 1 : length(time)
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
    
    engine = jtree_inf_engine(bnet{j});
    
    % set conditional nodes
    evidence = cell(1,k);
     
    evidence{index.mdo0}     = do0;
    evidence{index.mp0}      = p0;
    evidence{index.md2t(1)}  = d2t(1);
    evidence{index.md2t(2)}  = d2t(2);
    evidence{index.mf(1)}    = f1';
    evidence{index.mu(1)}    = u1';
    evidence{index.mf(3)}    = f3';
    evidence{index.mu(3)}    = u3';
    evidence{index.mdo(2)}   = do2';
    evidence{index.mp(2)}    = p2';
    
    [engine, ll] = enter_evidence(engine, evidence);
    
    known_do0(:,j) = do0;
    known_p0(:,j)  = p0;
    known_d2t(:,j) = d2t;
    known_f3(:,j)  = f3;
    known_u3(:,j)  = u3;
    
%     marg = marginal_nodes(engine, index.do(1));
%     footInertialHat(1:3,j) = reshape(marg.mu,3,1);    
%     marg = marginal_nodes(engine, index.p(1));
%     pMu(:,1) = reshape(marg.mu,3,1);
    
    marg = marginal_nodes(engine, index.p(2));
    footInertialHat(1:3,j) = reshape(marg.mu,3,1);
    footInertialSigma(1:3,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.do(2));
    footInertialHat(4:6,j) = reshape(marg.mu,3,1);
    footInertialSigma(4:6,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.f(1));
    legRightFTHat(1:3,j) = reshape(marg.mu,3,1);
    legRightFTSigma(1:3,j) = diag(reshape(marg.Sigma,3,3));
    
    marg = marginal_nodes(engine, index.u(1));
    legRightFTHat(4:6,j) = reshape(marg.mu,3,1);
    legRightFTSigma(4:6,j) = diag(reshape(marg.Sigma,3,3));

%     marg = marginal_nodes(engine, [index.f(2)]);
%     fMu(:,2) = reshape(marg.mu,3,1);
%     
%     marg = marginal_nodes(engine, [index.u(1)]);
%     uMu(:,1) = reshape(marg.mu,3,1);
%     
%     marg = marginal_nodes(engine, [index.c(1)]);
%     cMu(:,1) = reshape(marg.mu,3,1)
%     
%     marg = marginal_nodes(engine, [index.c(2)]);
%     cMu(:,2) = reshape(marg.mu,3,1)
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
% 
% figure
% subplot(221)
% plot(known_do0')
% subplot(222)
% plot(known_p0')
% subplot(223)
% plot(known_f3')
% subplot(224)
% plot(known_u3')
% 
% figure
% plot(known_d2t')

