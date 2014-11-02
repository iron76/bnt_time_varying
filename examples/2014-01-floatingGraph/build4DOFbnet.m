function [bnet, R] = build4DOFbnet(bnet, index, q, dq, o0)

defineConst;

c1 = cos(q(3));  c2 = cos(q(4));
s1 = sin(q(3));  s2 = sin(q(4));

% R23
R(1:3,1:3,3) = [s1 c1 0;
      -c1  s1 0;
      0    0 1];
% R34
R(1:3,1:3,4) = [c2 -s2 0;
      s2  c2 0;
      0    0 1];

R(1:3,1:3,5) = [1 0 0;
      0 1 0;
      0 0 1];  
% R01
R(1:3,1:3,1) = [1  0 0;
                0  0 1;
                0 -1 0];
% R12  
R(1:3,1:3,2) = [0 0 -1;
                0 1  0;
                1 0  0];

dt  = [dq(1), dq(2),  dq(3),  dq(4)];            

%Kinematic recursion
for i = 1:4
    Ri = R(:,:,i);
    if i==1
        o_prev = o0;
        index_do_prev = index.do0;
        index_p_prev  =  index.p0;
    else
        o_prev = o(:,i-1);        
        index_do_prev = index.do(i-1);
        index_p_prev  =  index.p(i-1);
    end
    
    if i == 3 || i == 4
        
        o(:,i)  = Ri'*(o_prev+dq(i)*z0);
        
        %% do(:,i) = Ri'*(do_prev+d2q(i)*z0+dq(i)*cross(o_prev,z0));
        b = Ri'*dq(i)*cross(o_prev,z0);
        W = commuteOnOrder2(Ri'*z0, Ri', index.d2t(i), index_do_prev);
        bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %%  p(:,i)  = Ri'*p_prev+cross(do(:,i),r0(:,i))+cross(o(:,i),cross(o(:,i),r0(:,i)));
        b = cross(o(:,i),cross(o(:,i),r0(:,i)));
        W = commuteOnOrder2(Ri', -vec_hat(r0(:,i)), index_p_prev, index.do(i));
        bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %% c(:,i)  = p(:,i) + cross(do(:,i), rc(:,i)) +cross(o(:,i),cross(o(:,i),rc(:,i)));
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        W = commuteOnOrder2(eye(d), -vec_hat(rc(:,i)), index.p(i), index.do(i));
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    else
        o(:,i)  = Ri'*o_prev;
        
        %% do(:,i) = Ri'*do_prev;
        b = zeros(size(o(:,i)));
        W = Ri';
        bnet.CPD{index.do(i)}  = gaussian_CPD(bnet, index.do(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);

        %% p(:,i)  = Ri'*(p_prev+d2t(i)*z0)+2*dt(i)*cross(o(:,i), Ri'*z0)+cross(do(:,i),ri0_i(:,i))+cross(o(:,i),cross(o(:,i),ri0_i(:,i)));
        b = 2*dt(i)*cross(o(:,i), Ri'*z0)+cross(o(:,i),cross(o(:,i),r0(:,i)));
        W = commuteOnOrder3(Ri', Ri'*z0, -vec_hat(r0(:,i)), index_p_prev, index.d2t(i), index.do(i));
        bnet.CPD{index.p(i)}  = gaussian_CPD(bnet, index.p(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
                
        %% c(:,i)  = p(:,i) + cross(do(:,i), riC_i(:,i)) +cross(o(:,i),cross(o(:,i),riC_i(:,i)));
        b = cross(o(:,i),cross(o(:,i),rc(:,i)));
        W = commuteOnOrder2(eye(d), -vec_hat(rc(:,i)), index.p(i), index.do(i));
        bnet.CPD{index.c(i)}  = gaussian_CPD(bnet, index.c(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
    end
end

Ree = -1.*(Rb0*R(:,:,1)*R(:,:,2)*R(:,:,3)*R(:,:,4)*R(:,:,5))';
% f(5) = Ree * fee;
% u(5) = Ree * uee;
bnet.CPD{index.f(5)}  = gaussian_CPD(bnet, index.f(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', Ree);
bnet.CPD{index.u(5)}  = gaussian_CPD(bnet, index.u(5),  'mean', b, 'cov', sModel.*eye(d), 'weights', Ree);

%Dynamic recursion
for i = 4:-1:1
    Ri = R(:,:,i+1);
    Ii = I(:,:,i);
    if i == 3
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i) + Rbb*fbb;
        %  fbbNE   = Rbb * fbb;
        %  ubbNE   = Rbb * ubb;
        %  Rbb     = -1.*(Rb0*R(:,:,1)*R(:,:,2)*R(:,:,3))';
        
        b = zeros(d,1);
        Rbb = -1.*(Rb0*R(:,:,1)*R(:,:,2)*R(:,:,3))';
        W = commuteOnOrder3(Ri, eye(d).*m(i,1), Rbb, index.f(i+1), index.c(i), index.fbb);
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i)) + cross(-rc(:,i)-r0(:,i), Rbb*fbb) + Rbb*ubb;
        b = cross(o(:,i), Ii*o(:,i));
        Wt(:, :, 1) = Ii;
        Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
        Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
        Wt(:, :, 4) = Ri;
        Wt(:, :, 5) = vec_hat(rc(:,i)+r0(:,i))*(Rb0*R(:,:,1)*R(:,:,2)*R(:,:,3))';
        Wt(:, :, 6) = Rbb;
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1), index.fbb, index.ubb]);
        W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4)) Wt(:,:,ind(5)) Wt(:,:,ind(6))];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                
    else
        %% f(:,i)  = Ri*f(:,i+1) + m(i,1).*c(:,i);
        b = zeros(d,1);
        W = commuteOnOrder2(Ri, eye(d).*m(i,1), index.f(i+1), index.c(i));
        bnet.CPD{index.f(i)}  = gaussian_CPD(bnet, index.f(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);
        
        %% u(:,i)  = -cross(f(:,i), rc(:,i)+r0(:,i)) + Ri*u(:,i+1)+ ...
        %    cross(Ri*f(:,i+1), rc(:,i)) + Ii*do(:,i) + cross(o(:,i), Ii*o(:,i))
        b = cross(o(:,i), Ii*o(:,i));
        Wt(:, :, 1) = Ii;
        Wt(:, :, 2) = vec_hat(rc(:,i)+r0(:,i));
        Wt(:, :, 3) = -vec_hat(rc(:,i))*Ri;
        Wt(:, :, 4) = Ri;
        [~,ind] = sort([index.do(i), index.f(i), index.f(i+1),  index.u(i+1)]);
        W = [Wt(:,:,ind(1)) Wt(:,:,ind(2)) Wt(:,:,ind(3)) Wt(:,:,ind(4))];
        
        bnet.CPD{index.u(i)}  = gaussian_CPD(bnet, index.u(i),  'mean', b, 'cov', sModel.*eye(d), 'weights', W);                

    end
end

