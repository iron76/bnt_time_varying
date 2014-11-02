function bnet = build2DOFbnet(bnet, index, t, dt, d2t, o0)

defineConst2

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
        W = commuteOnOrder(Ri'*z0, Ri', -eye(3), index.d2t(i), index.do0, index.do(i));
        bnet.CPD{index.edo(i)}  = gaussian_CPD(bnet, index.edo(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        W = commuteOnOrder(Ri', -vec_hat(r0(:,i)), -eye(3), index.p0, index.do(i), index.p(i));
        bnet.CPD{index.ep(i)}  = gaussian_CPD(bnet, index.ep(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        W = commuteOnOrder(eye(d), -vec_hat(rc(:,i)), -eye(3), index.p(i), index.do(i), index.c(i));
        bnet.CPD{index.ec(i)}  = gaussian_CPD(bnet, index.ec(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
    else
        o(:,i)  = Ri'*(o(:,i-1)+dt(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2t(i)*z0+dt(i)*cross(o_prev,z0));
        b = Ri'*dt(i)*cross(o(:,i-1),z0);
        W = commuteOnOrder(Ri'*z0, Ri', -eye(3), index.d2t(i), index.do(i-1), index.do(i));
        bnet.CPD{index.edo(i)}  = gaussian_CPD(bnet, index.edo(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        W = commuteOnOrder(Ri', -vec_hat(r0(:,i)), -eye(3), index.p(i-1), index.do(i), index.p(i));
        bnet.CPD{index.ep(i)}  = gaussian_CPD(bnet, index.ep(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%
        % c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        W = commuteOnOrder(eye(d), -vec_hat(rc(:,i)), -eye(3), index.p(i), index.do(i), index.c(i));
        bnet.CPD{index.ec(i)}  = gaussian_CPD(bnet, index.ec(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
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
    W = commuteOnOrder(Ri, eye(d).*m(i,1), -eye(3), index.f(i+1), index.c(i), index.f(i));
    bnet.CPD{index.ef(i)}  = gaussian_CPD(bnet, index.ef(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    
    %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
    %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i));
    b = cross(o(:,i), Ii*o(:,i));
    Wt(:, :, 1) = Ii;
    Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
    Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
    Wt(:, :, 4) = Ri;
    Wt(:, :, 5) = -eye(3);
    [y,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1), index.u(i)]);
    W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))  Wt(:,:,ind(5))];
    
    bnet.CPD{index.eu(i)}  = gaussian_CPD(bnet, index.eu(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
end
