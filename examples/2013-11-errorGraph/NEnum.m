clear all
close all

defineConst2;
dtime = 0.40;                %time discretization [s]
[footInertial, tm] = readDataDumper('../../data/dump2/inertial/data.log');
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



for j = 1 : length(time)
    f3 = footRightFT(j,1:3);
    u3 = footRightFT(j,4:6);
    
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
    
    f(:,3) = f3;
    u(:,3) = u3;
    
    %Useful definitions
    c1 = cos(t(1));  c2 = cos(t(2));
    s1 = sin(t(1));  s2 = sin(t(2));
    
    R_10 = [c1 -s1 0;
        s1  c1 0;
        0    0 1];
    
    R_21 = [c2 -s2 0;
        s2  c2 0;
        0    0 1];
    
    R_32 = [1 0 0;
        0 1 0;
        0 0 1];
    
    dt  = [dt(1),  dt(2)];
    d2t = [d2t(1), d2t(2)];
    
    %Kinematic recursion
    for i = 1:2
        if i==1
            Ri=R_10;
        else
            Ri=R_21;
        end
        if i==1
            o_prev = o0;
            p_prev = p0;
            do_prev = do0;
        else
            p_prev = p(:,i-1);
            o_prev = o(:,i-1);
            do_prev = do(:,i-1);
        end
        o(:,i)  = Ri'*(o_prev+dt(i)*z0);
        do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
    end
    
    %Dynamic recursion
    for i = 2:-1:1
        if i==2
            Ri=R_32;
            Ii = I2;
        else
            Ri=R_21;
            Ii = I1;
        end
        f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i);
        u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
            cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i));
    end
    footInertialHat(1:3,j) = p(:,2);
    footInertialHat(4:6,j) = do(:,2);
    legRightFTHat(1:3,j) = f(:,1);
    legRightFTHat(4:6,j) = u(:,1);
end

%inertial sensor prediction
subplot(221)
plot(time, footInertialHat(1:3,:)')
title('estimated')
est_ax = axis;
subplot(222)
plot(time, -footInertial(:,4:6))     %d2p_computed = -d2p_meas
title('measured')
axis(est_ax);
subplot(223)
plot(time, footInertialHat(4:6,:)')
est_ax = axis;
subplot(224)
plot(time(1:end-1), diff(footInertial(:,7:9)))
axis(est_ax);

figure
%FT sensor prediction
subplot(221)
plot(time, legRightFTHat(1:3,:)')
title('estimated')
est_ax = axis;
subplot(222)
plot(time, legRightFT(:,1:3))
title('measured')
axis(est_ax);
subplot(223)
plot(time, legRightFTHat(4:6,:)')
est_ax = axis;
subplot(224)
plot(time, legRightFT(:,4:6))
axis(est_ax);

