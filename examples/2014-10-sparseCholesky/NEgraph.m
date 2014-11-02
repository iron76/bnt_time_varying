clear all
close all
defineConst;

dtime = 0.40;                %time discretization [s]
[footInertial, tm] = readDataDumper('../../data/dump2/inertial/data.log');
[footInertial(:,7), dfootInertial(:,1)] = sgolayDiff(footInertial(:,7), tm);
[footInertial(:,8), dfootInertial(:,2)] = sgolayDiff(footInertial(:,8), tm);
[footInertial(:,9), dfootInertial(:,3)] = sgolayDiff(footInertial(:,9), tm);
footInertial(:,7:9)  = pi/180.*footInertial(:,7:9);
dfootInertial(:,1:3) = pi/180.*(Adj_IMU(1:3, 1:3) * dfootInertial(:,1:3)')';
time = tm(300):dtime:tm(end-300);    %time samples [s]
footInertial  = interp1(tm,  footInertial, time, 'spline', 'extrap');
dfootInertial = interp1(tm, dfootInertial, time, 'spline', 'extrap');
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

[bnet, index, new_ns] = init4DOFnet(n, k, d);
for j = 1 : length(time)
    %Useful definitions

    f5 = footRightFT(j,1:3);
    u5 = footRightFT(j,4:6);

    f3 = legRightFT(j,1:3);
    u3 = legRightFT(j,4:6);
    
    o4  =  footInertial(j,7:9);
    p4  = -footInertial(j,4:6);
        
    do4 =  dfootInertial(j,1:3);
    
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
    
    fbb = [0 0 0]';
    ubb = [0 0 0]';
    
    y = [...
        p4';                    o4'; ...
        f3';                    u3'; ...        
        f5';                    u5'; ...
        d2t(3);                 d2t(4); ... 
        fbb;                    ubb;];
    
    
    [h, Sh, R, avgTime(j)]  = build4DOFestim(index, y, t, dt, o0, do0, p0);
    
    if j == 1
        [h, Sh, R, avgTime2(j), S] = build4DOFchol(index, y, t, dt, o0, do0, p0, j);
    else
        [h, Sh, R, avgTime2(j)] = build4DOFchol(index, y, t, dt, o0, do0, p0, j, S);
    end
    
    
    footInertialHat(1:3,j) = reshape(h{index.p(4)},3,1);
    footInertialSigma(1:3,j) = diag(reshape(Sh{index.p(4)},3,3));
    
%     footInertialHat(1:3,j) 
%     Exy(1:3)
    
    footInertialHat(4:6,j) = reshape(h{index.o(4)},3,1);
    footInertialSigma(4:6,j) = diag(reshape(Sh{index.o(4)},3,3));
    
    legRightFTHat(1:3,j) = R(1:3,1:3,3)'*reshape(h{index.f(2)},3,1);
    legRightFTSigma(1:3,j) = diag(reshape(Sh{index.f(2)},3,3));
    
    legRightFTHat(4:6,j) = R(1:3,1:3,3)'*reshape(h{index.u(2)},3,1);
    legRightFTSigma(4:6,j) = diag(reshape(Sh{index.u(2)},3,3));

    legRightd2tHat(1:4,j) = [h{index.d2t(1)}; h{index.d2t(2)}; h{index.d2t(3)}; h{index.d2t(4)}];
    legRightd2tSigma(1:4,j) = [Sh{index.d2t(1)}; Sh{index.d2t(2)}; Sh{index.d2t(3)}; Sh{index.d2t(4)}];
            
    fprintf(1, 'total time is: %% %.1f in %.5f sec VS %.5f sec\r ', j/length(time)*100, avgTime(j), avgTime2(j))
end

close all

%inertial sensor prediction
subplot(221)
plotWithSigma(time, footInertialHat(1:3,:)', footInertialSigma(1:3,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plotWithSigma(time, -footInertial(:,4:6), repmat(diag(Sp4)', length(footInertial(:,4)), 1), '')     %d2p_computed = -d2p_meas
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, footInertialHat(4:6,:)', footInertialSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plotWithSigma(time(1:end), footInertial(:,7:9), repmat(diag(So4)', length(footInertial(:,7)), 1), '')
axis(est_ax);

figure
%FT sensor prediction
subplot(221)
plotWithSigma(time, legRightFTHat(1:3,:)', legRightFTSigma(1:3,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plotWithSigma(time, legRightFT(:,1:3), repmat(diag(Sf2)', length(legRightFT(:,1)), 1), '')
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, legRightFTHat(4:6,:)', legRightFTSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plotWithSigma(time, legRightFT(:,4:6), repmat(diag(Su2)', length(legRightFT(:,4)), 1), '')
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
