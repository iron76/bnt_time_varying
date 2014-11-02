clear all
close all

defineConstState;

dtime = 0.05;                %time discretization [s]
[footInertial, tm] = readDataDumper('../../data/dump2/inertial/data.log');
% footInertial(:, 7) = footInertial(:, 7) - mean(footInertial(1:1000, 7));
% footInertial(:, 8) = footInertial(:, 8) - mean(footInertial(1:1000, 8));
% footInertial(:, 9) = footInertial(:, 9) - mean(footInertial(1:1000, 9));
[footInertial(:,7), dfootInertial(:,1)] = sgolayDiff(footInertial(:,7), tm);
[footInertial(:,8), dfootInertial(:,2)] = sgolayDiff(footInertial(:,8), tm);
[footInertial(:,9), dfootInertial(:,3)] = sgolayDiff(footInertial(:,9), tm);
footInertial(:,7:9)  = pi/180.*footInertial(:,7:9);
dfootInertial(:,1:3) = pi/180.*(Adj_IMU(1:3, 1:3) * dfootInertial(:,1:3)')';
time = tm(200):dtime:tm(end-200);    %time samples [s]
time = time(1,140:end);
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


bnet = cell(length(time), 1);
[bnet0, index, labels] = init4DOFnetState(n, k, d);

x = [0 0 0 0, 0 0 0 0]';
P = diag([1e4 1e4 1e-4 1e-2, 1e2 1e2 1e-2 1e-2]);
engine = [];

%dfootInertial(:,1:3) = 0;
%footInertial(:,7:9) = 0;


for j = 1 : length(time)
    
    %Useful definitions
    f5 = footRightFT(j,1:3);
    u5 = footRightFT(j,4:6);

    f2 = legRightFT(j,1:3);
    u2 = legRightFT(j,4:6);
    
    do4 = dfootInertial(j,1:3);
    o4  =  footInertial(j,7:9);
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
    d2q = [0 0 legRight_d2q(j,[1 4])*pi/180];
    
    fbb = [0 0 0]';
    ubb = [0 0 0]';
    
    bnet{j}      = bnet0;
    
    y = [d2q(4); f2'; u2'; f5'; u5'; o4'; p4'; [0 0 0]'; [0 0 0]'];
    if j == 1
        [h, Sh, indexH, R] = computeHidden(t, dt, y, o0, p0, do0);
        x(3) = t(3);  x(7) = dt(3); 
    end
    z = [p0;               do0; ...
        h{indexH.p(1)};   h{indexH.p(2)};  h{indexH.p(3)};  p4';  ...
        h{indexH.do(1)};  h{indexH.do(2)}; h{indexH.do(3)}; do4'; ...
        h{indexH.f(1)};   R(1:3,1:3,3)*f2';h{indexH.f(3)};  h{indexH.f(4)};   f5'; ...
        h{indexH.u(1)};   R(1:3,1:3,3)*u2';h{indexH.u(3)};  h{indexH.u(4)};   u5'; ...
        h{indexH.d2q(1)}; h{indexH.d2q(2)};h{indexH.d2q(3)};          d2q(4); ...
        h{indexH.c(1)};   h{indexH.c(2)};  h{indexH.c(3)};  h{indexH.c(4)};  ...
        fbb;              ubb;             h{indexH.fee};   h{indexH.uee};    ...
        h{indexH.o(1)};   h{indexH.o(2)};  h{indexH.o(3)};  o4'];
    
    x(4) = t(4);    %x(3) = t(3);    
    x(8) = dt(4);   %x(7) = dt(3);
    % x = [0 0 t(3) t(4), 0 0 dt(3) dt(4)]';
    % P = diag([1e2 1e2 1e2 1e-3, 1e2 1e2 1e2 1e-3]);
    [h, Sh, indexH, R, engine]    = build4DOFbnetState(bnet{j}, engine, index, x, P, z, o0, p0, do0);
         
    footInertialHat(1:3,j) = reshape(h{indexH.p(4)},3,1);
    footInertialSigma(1:3,j) = diag(reshape(Sh{indexH.p(4)},3,3));
    
    footInertialHat(4:6,j) = reshape(h{indexH.o(4)},3,1);
    footInertialSigma(4:6,j) = diag(reshape(Sh{indexH.o(4)},3,3));
    
    legRightFTHat(1:3,j) = R(1:3,1:3,3)'*reshape(h{indexH.f(2)},3,1);
    legRightFTSigma(1:3,j) = diag(reshape(Sh{indexH.f(2)},3,3));
    
    legRightFTHat(4:6,j) = R(1:3,1:3,3)'*reshape(h{indexH.u(2)},3,1);
    legRightFTSigma(4:6,j) = diag(reshape(Sh{indexH.u(2)},3,3));

    legRightd2qHat(1:4,j) = [h{indexH.d2q(1)}; h{indexH.d2q(2)}; h{indexH.d2q(3)}; h{indexH.d2q(4)}];
    legRightd2qSigma(1:4,j) = [Sh{indexH.d2q(1)}; Sh{indexH.d2q(2)}; Sh{indexH.d2q(3)}; Sh{indexH.d2q(4)}];

    
    newo1(:,j) = h{indexH.o(1)};
    newo2(:,j) = h{indexH.o(2)};

    newo3(:,j) = h{indexH.o(3)};
    newo4(:,j) = h{indexH.o(4)};
    
    xHat = h{indexH.x};
    PHat = Sh{indexH.x};
   
    fprintf(1, 'total time is: %% %.1f\r', j/length(time)*100)
    
    A    = [eye(n)       eye(n)*dtime; zeros(n) eye(n)];
    Q    = [eye(n)/10*2e-2   zeros(n); zeros(n) eye(n)*dtime/0.2];
    Q(7,7) = 5e-3;   Q(8,8) = 1e-5;
    Q(3,3) = 1e-5;
    PHat = A*PHat*A' +  Q;
    % xHat = A*xHat;
    %(Sh{indexH.o(3)}(3,3)-PHat(7,7))/PHat(7,7)
    
    xHat(3) = xHat(3) + xHat(7)*dtime;
    xHat(4) = xHat(4) + xHat(8)*dtime;
    legRightdqHat(1:4,j)   = [xHat(5)   xHat(6)   xHat(7)   xHat(8)];
    legRightdqSigma(1:4,j) = [PHat(5,5) PHat(6,6) PHat(7,7) PHat(8,8)];

    legRightqHat(1:4,j)    = [xHat(1)   xHat(2)   xHat(3)   xHat(4)];
    legRightqSigma(1:4,j)  = [PHat(1,1) PHat(2,2) PHat(3,3) PHat(4,4)];    

    x = xHat;
    P = PHat;    
end

close all
D = bnet0.dag;
%remove excessive connection to x
D(index.x,:) = 0;
draw_graph(D, labels);

%inertial sensor prediction
figure
subplot(221)
plotWithSigma(time, footInertialHat(1:3,:)', footInertialSigma(1:3,:)', '')
title('estimated lin acc and ang vel')
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
title('FT')
%FT sensor prediction
subplot(221)
plotWithSigma(time, legRightFTHat(1:3,:)', legRightFTSigma(1:3,:)', '')
title('estimated Force and Torques')
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
title('d2q')
%d2q sensor prediction
subplot(321)
plotWithSigma(time, legRightd2qHat(1:2,:)', legRightd2qSigma(1:2,:)', '')
title('estimated base and joint accelerations')
est_ax = axis;
subplot(322)
plotWithSigma(time, 0.*legRight_d2q(:,[1 4])*pi/180, repmat([Sd2q1 Sd2q2], length(legRight_d2q(:,[1 4])), 1), '')
axis(est_ax);
title('measured')
subplot(323)
plotWithSigma(time, legRightd2qHat(3,:)', legRightd2qSigma(3,:)', '')
est_ax = axis;
subplot(324)
plotWithSigma(time, legRight_d2q(:,[1])*pi/180, repmat([Sd2q3], length(legRight_d2q(:,[1 4])), 1), '')
axis(est_ax);
subplot(325)
plotWithSigma(time, legRightd2qHat(4,:)', legRightd2qSigma(4,:)', '')
est_ax = axis;
subplot(326)
plotWithSigma(time, legRight_d2q(:,[4])*pi/180, repmat([Sd2q4], length(legRight_d2q(:,[1 4])), 1), '')
axis(est_ax);


figure
title('dq')
%dq sensor prediction
subplot(321)
plotWithSigma(time, legRightdqHat(1:2,:)', legRightdqSigma(1:2,:)', '')
title('estimated base and joint velocities')
est_ax = axis;
subplot(322)
plotWithSigma(time, 0.*legRight_dq(:,[1 4])*pi/180, repmat([P(5,5) P(6,6)], length(legRight_dq(:,[1 4])), 1), '')
axis(est_ax);
title('measured')
subplot(323)
plotWithSigma(time, legRightdqHat(3,:)', legRightdqSigma(3,:)', '')
est_ax = axis;
subplot(324)
plotWithSigma(time, legRight_dq(:,[1])*pi/180, repmat(Sdq4, length(legRight_dq(:,[1 4])), 1), '')
axis(est_ax);
subplot(325)
plotWithSigma(time, legRightdqHat(4,:)', legRightdqSigma(4,:)', '')
est_ax = axis;
subplot(326)
plotWithSigma(time, legRight_dq(:,[4])*pi/180, repmat(Sdq4, length(legRight_dq(:,[1 4])), 1), '')
axis(est_ax);


figure
title('q')
%dq sensor prediction
subplot(321)
plotWithSigma(time, legRightqHat(1:2,:)', legRightqSigma(1:2,:)', '')
title('estimated base and joint positions')
est_ax = axis;
subplot(322)
plotWithSigma(time, 0.*legRight_q(:,[1 4])*pi/180, repmat([P(1,1) P(2,2)], length(legRight_q(:,[1 4])), 1), '')
axis(est_ax);
title('measured')
axis(est_ax);
subplot(323)
plotWithSigma(time, legRightqHat(3,:)', legRightqSigma(3,:)', '')
est_ax = axis;
subplot(324)
plotWithSigma(time, legRight_q(:,[1])*pi/180, repmat(Sq4, length(legRight_q(:,[1 4])), 1), '')
axis(est_ax);
axis(est_ax);
subplot(325)
plotWithSigma(time, legRightqHat(4,:)', legRightqSigma(4,:)', '')
est_ax = axis;
subplot(326)
plotWithSigma(time, legRight_q(:,[4])*pi/180, repmat(Sq4, length(legRight_q(:,[1 4])), 1), '')
axis(est_ax);
