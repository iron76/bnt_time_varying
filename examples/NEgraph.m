clear all
close all

defineConst;

dtime = 0.40;                %time discretization [s]
[footInertial, tm] = readDataDumper('../data/dump2/inertial/data.log');
time = tm(200):dtime:tm(end-200);    %time samples [s]
footInertial = interp1(tm, footInertial, time, 'spline', 'extrap');
footInertial(:,4:9) = (Adj_IMU*footInertial(:,4:9)')';

[footRightFT, tm] = readDataDumper('../data/dump2/rightFootFT/data.log');
footRightFT = interp1(tm, footRightFT, time, 'spline', 'extrap');
footRightFT = (Adj_FT*footRightFT')';

[footRightSkin, tm] = readDataDumper('../data/dump2/rightFootSkin/data.log');
footRightSkin = interp1(tm, footRightSkin, time, 'spline', 'extrap');

[legRightFT, tm] = readDataDumper('../data/dump2/rightLegFT/data.log');
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

[legRight_q, tm] = readDataDumper('../data/dump2/rightLegPos/data.log');
[tm, I, J] = unique(tm);
legRight_q = legRight_q(I, :);
legRight_q = interp1(tm, legRight_q, time, 'spline', 'extrap');
legRight_q(:,1) = -legRight_q(:,1)-q0_offset(1);
legRight_q(:,4) = -legRight_q(:,4)-q0_offset(4);

[legRight_dq, tm] = readDataDumper('../data/dump2/rightLegVel/data.log');
[tm, I, J] = unique(tm);
legRight_dq = legRight_dq(I, :);
legRight_dq = interp1(tm, legRight_dq, time, 'spline', 'extrap');

[legRight_d2q, tm] = readDataDumper('../data/dump2/rightLegAcc/data.log');
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

    f1 = legRightFT(j,1:3);
    u1 = legRightFT(j,4:6);
    
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
    
    [bnet{j}, index] = init2DOFnet(n, k, d);
    
    bnet{j}.CPD{index.do0}    = gaussian_CPD(bnet{j}, index.do0,    'mean', do0,    'cov', Sdo0);
    bnet{j}.CPD{index.p0}     = gaussian_CPD(bnet{j}, index.p0,     'mean', p0,     'cov', Sp0);
    bnet{j}.CPD{index.d2t(1)} = gaussian_CPD(bnet{j}, index.d2t(1), 'mean', d2t(1), 'cov', Sd2t1);
    bnet{j}.CPD{index.d2t(2)} = gaussian_CPD(bnet{j}, index.d2t(2), 'mean', d2t(2), 'cov', Sd2t2);
    bnet{j}.CPD{index.f(3)}   = gaussian_CPD(bnet{j}, index.f(3),   'mean', f3,     'cov', Sf3);
    bnet{j}.CPD{index.u(3)}   = gaussian_CPD(bnet{j}, index.u(3),   'mean', u3,     'cov', Su3);
    
    %Useful definitions
    c1 = cos(t(1));  c2 = cos(t(2));
    s1 = sin(t(1));  s2 = sin(t(2));
    
    R1 = [c1 -s1 0;
        s1  c1 0;
        0    0 1];
    
    R2 = [c2 -s2 0;
        s2  c2 0;
        0    0 1];
    
    R3 = [1 0 0;
        0 1 0;
        0 0 1];
    
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
            if (index.d2t(i) < index.do0)
                W = [Ri'*z0 Ri'];
            else
                W = [Ri' Ri'*z0];
            end
            bnet{j}.CPD{index.do(i)}  = gaussian_CPD(bnet{j}, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            
            %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
            b = cross(o(:,i),cross(o(:,i),r0(:,i)));
            if (index.p0 < index.do(i))
                W = [Ri' -vec_hat(r0(:,i))];
            else
                W = [-vec_hat(r0(:,i)) Ri'];
            end
            bnet{j}.CPD{index.p(i)}  = gaussian_CPD(bnet{j}, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            
            %%
            % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
            b = cross(o(:,i),cross(o(:,i),rc(:,i)));
            if (index.p(i) < index.do(i))
                W = [eye(d) -vec_hat(rc(:,i))];
            else
                W = [-vec_hat(rc(:,i)) eye(d)];
            end
            bnet{j}.CPD{index.c(i)}  = gaussian_CPD(bnet{j}, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            
        else
            o(:,i)  = Ri'*(o(:,i-1)+dt(i)*z0);
            
            %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
            b = Ri'*dt(i)*cross(o(:,i-1),z0);
            if (index.d2t(i) < index.do(i-1))
                W = [Ri'*z0 Ri'];
            else
                W = [Ri' Ri'*z0];
            end
            bnet{j}.CPD{index.do(i)}  = gaussian_CPD(bnet{j}, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            
            %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
            b = cross(o(:,i),cross(o(:,i),r0(:,i)));
            if (index.p(i-1) < index.do(i))
                W = [Ri' -vec_hat(r0(:,i))];
            else
                W = [-vec_hat(r0(:,i)) Ri'];
            end
            bnet{j}.CPD{index.p(i)}  = gaussian_CPD(bnet{j}, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
            
            %%
            % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
            b = cross(o(:,i),cross(o(:,i),rc(:,i)));
            if (index.p(i) < index.do(i))
                W = [eye(d) -vec_hat(rc(:,i))];
            else
                W = [-vec_hat(rc(:,i)) eye(d)];
            end
            bnet{j}.CPD{index.c(i)}  = gaussian_CPD(bnet{j}, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
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
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i);
        b = zeros(d,1);
        if (index.f(i+1) < index.c(i))
            W = [Ri eye(d).*m(i,1)];
        else
            W = [eye(d).*m(i,1) Ri];
        end
        bnet{j}.CPD{index.f(i)}  = gaussian_CPD(bnet{j}, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i));
        b = cross(o(:,i), Ii*o(:,i));
        
        Wt(:, :, 1) = Ii;
        Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
        Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
        Wt(:, :, 4) = Ri;
        [y,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1)]);
        W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))];
        
        bnet{j}.CPD{index.u(i)}  = gaussian_CPD(bnet{j}, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    end
    
    engine = jtree_inf_engine(bnet{j});
    
    % set conditional nodes
    evidence = cell(1,k);
%     evidence{index.do0   }  = do0;
%     evidence{index.p0    }  = p0;
%     evidence{index.d2t(1)}  = d2t(1);
%     evidence{index.d2t(2)}  = d2t(2);
%     evidence{index.f(3)  }  = f3;
%     evidence{index.u(3)  }  = u3;

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
plot(time, -footInertial(:,4:6))     %d2p_computed = -d2p_meas
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, footInertialHat(4:6,:)', footInertialSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plot(time(1:end-1), diff(footInertial(:,7:9)))
axis(est_ax);

figure
%FT sensor prediction
subplot(221)
plotWithSigma(time, legRightFTHat(1:3,:)', legRightFTSigma(1:3,:)', '')
title('estimated')
est_ax = axis;
subplot(222)
plot(time, legRightFT(:,1:3))
title('measured')
axis(est_ax);
subplot(223)
plotWithSigma(time, legRightFTHat(4:6,:)', legRightFTSigma(4:6,:)', '')
est_ax = axis;
subplot(224)
plot(time, legRightFT(:,4:6))
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

