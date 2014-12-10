function  [d, Sd, net] = nID( model, q, qd, qdd, f_ext, y, net)

% ID  Inverse Dynamics via Recursive Newton-Euler Algorithm
% ID(model,q,qd,qdd,f_ext) calculates the inverse dynamics of a kinematic
% tree via the recursive Newton-Euler algorithm.  q, qd and qdd are vectors
% of joint position, velocity and acceleration variables; and the return
% value is a vector of joint force variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

% n  = model.NB;
% Size of the D matrix [dm, dn]
% dm = 18*model.NB + n
% dn = 24 * model.NB + 2*n

global Sm Su

NB  = model.NB;

[ny, ~]   = size(y.values);

a_grav = get_gravity(model);

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  [~, jn{i}] = size(S{i});
  for j = 1:model.NB
    Dc{i,j} = zeros(18+jn{i}, 24+2*jn{i});
  end
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    % a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd(i) + crm(v{i})*vJ;
  end
end

if nargin >= 5
  [~, Xa] = apply_external_forces( model.parent, Xup, cell(NB), cell(NB) );
end


for i = 1:model.NB
  if model.parent(i) == 0
    % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
    ai   = net.nodes.index{i}(1);
    nai  = net.nodes.sizes{i,1}(1);
    d2qi = net.nodes.index{i}(6);
    pars  = cell2mat(net.bnet.parents(ai));
    Wa    = cell(1, length(pars));
    for j = 1 : length(pars)
      if pars(j) == d2qi
        Wa{1, j} = S{i};
      end
    end
    W = cell2mat(Wa);
    net.bnet.CPD{ai} = gaussian_CPD(net.bnet, ai, 'mean', Xup{i}*(-a_grav), 'cov', Sm.a{i}, 'weights', W);      
  else
    % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
    vJ = S{i}*qd(i);
    
    % a{i} = Xup{i}*a{model.parent(i)} + ... 
    ai      = net.nodes.index{i}(1);
    d2qi    = net.nodes.index{i}(6);
    nai     = net.nodes.sizes{i,1}(1);
    aj      = net.nodes.index{model.parent(i)}(1);
    
    pars  = cell2mat(net.bnet.parents(ai));
    Wa    = cell(1, length(pars));
    for j = 1 : length(pars)
      if pars(j) == aj
        Wa{1, j} = Xup{i};
      elseif pars(j) == d2qi
        Wa{1, j} = S{i};
      end
    end
    W = cell2mat(Wa);
    net.bnet.CPD{ai} = gaussian_CPD(net.bnet, ai, 'mean', crm(v{i})*vJ, 'cov', Sm.a{i}, 'weights', W);
  end
  % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
  fBi  = net.nodes.index{i}(2);
  nfBi = net.nodes.sizes{i,1}(2);
  pars  = cell2mat(net.bnet.parents(fBi));
  Wa    = cell(1, length(pars));
  for j = 1 : length(pars)
    if pars(j) == ai
      Wa{1, j} = model.I{i};
    end
  end
  W = cell2mat(Wa);
  net.bnet.CPD{fBi} = gaussian_CPD(net.bnet, fBi, 'mean', crf(v{i})*model.I{i}*v{i}, 'cov', Sm.fB{i}, 'weights', W);    
end

for i = model.NB:-1:1
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
  fBi = net.nodes.index{i}(2);
  fi  = net.nodes.index{i}(3);
  fxi = net.nodes.index{i}(5);
  nfi = net.nodes.sizes{i,1}(3);
  pars  = cell2mat(net.bnet.parents(fi));
  Wa    = cell(1, length(pars));
  for j = 1 : length(pars)
    if pars(j) == fBi
      Wa{1, j} = eye(nfi);
    elseif pars(j) == fxi
      Wa{1, j} = -inv(Xa{i}');
    else
      Wa{1, j} = Xup{net.link(pars(j))}';
    end
  end
  W = cell2mat(Wa);
  
  net.bnet.CPD{fi} = gaussian_CPD(net.bnet, fi, 'mean', zeros(nfi,1), 'cov', Sm.f{i}, 'weights', W);    

  % tau(i,1) = S{i}' * f{i};
  taui  = net.nodes.index{i}(4);
  ntaui = net.nodes.sizes{i,1}(4);
  pars  = cell2mat(net.bnet.parents(taui));
  Wa    = cell(1, length(pars));
  for j = 1 : length(pars)
    if pars(j) == fi
      Wa{1, j} = S{i}';
    end
  end
  W = cell2mat(Wa);
  net.bnet.CPD{taui} = gaussian_CPD(net.bnet, taui, 'mean', zeros(ntaui, 1), 'cov', Sm.tau{i}, 'weights', W);    
  
  % fxi
  fxi  = net.nodes.index{i}(5);
  nfxi = net.nodes.sizes{i,1}(5);
  net.bnet.CPD{fxi} = gaussian_CPD(net.bnet, fxi, 'mean', zeros(nfxi, 1), 'cov', Su.fx{i});  
  
  % d2qi
  d2qi  = net.nodes.index{i}(6);
  nd2qi = net.nodes.sizes{i,1}(6);
  net.bnet.CPD{d2qi} = gaussian_CPD(net.bnet, d2qi, 'mean', zeros(nd2qi, 1), 'cov', Su.d2q{i});  
  
end

for i = 1 : ny
  yi  = net.nodes.index{NB*6 + i};
  nyi = net.nodes.sizes{NB*6 + i};
  net.bnet.CPD{yi} = gaussian_CPD(net.bnet, yi, 'mean', zeros(nyi, 1), 'cov', y.S{i,1}, 'weights', y.W{i,1});    
end

a   = zeros(NB,6);
fB  = zeros(NB,6);
f   = zeros(NB,6);
tau = zeros(NB,1);

Sa   = cell(NB,1);
SfB  = cell(NB,1);
Sf   = cell(NB,1);
Stau = cell(NB,1);
Sfx  = cell(NB,1);
Sd2q = cell(NB,1);

if isempty(net.engine)  
  net.engine = jtree_inf_engine(net.bnet);
else
  % disp('Using previous path on jtree!')
  net.engine = bnet_to_engine(net.bnet, net.engine);
end


evidence  = cell(1,6*NB+ny);
for i = 1 : ny
  evidence{net.nodes.index{NB*6 + i}} = y.values{i};
end
net.engine = enter_evidence(net.engine, evidence);
for i = NB:-1:1
  tmp        = marginal_nodes(net.engine, net.nodes.index{i}(1));
  a(i,:)   = tmp.mu';
  Sa{i}    = tmp.Sigma;
  tmp      = marginal_nodes(net.engine, net.nodes.index{i}(2));
  fB(i,:)  = tmp.mu';
  SfB{i}   = tmp.Sigma;
  tmp      = marginal_nodes(net.engine, net.nodes.index{i}(3));
  f(i,:)   = tmp.mu';
  Sf{i}    = tmp.Sigma;
  tmp      = marginal_nodes(net.engine, net.nodes.index{i}(4));
  tau(i,1) = tmp.mu;
  Stau{i}  = tmp.Sigma;
  tmp      = marginal_nodes(net.engine, net.nodes.index{i}(5));
  fx(i,:)  = tmp.mu';
  Sfx{i}   = tmp.Sigma;
  tmp      = marginal_nodes(net.engine, net.nodes.index{i}(6));
  d2q(i,1) = tmp.mu';
  Sd2q{i}  = tmp.Sigma;
  
end

for i = 1 : NB
  Sd((1:6)+(i-1)*19, (1:6)+(i-1)*19) = Sa{i};
  Sd((7:12)+(i-1)*19, (7:12)+(i-1)*19) = SfB{i};
  Sd((13:18)+(i-1)*19, (13:18)+(i-1)*19) = Sf{i};
  Sd((19:19)+(i-1)*19, (19:19)+(i-1)*19) = Stau{i};

  d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end

for i = 1 : NB
  Sd((1:6)+(i-1)*7+19*NB, (1:6)+(i-1)*7+19*NB) = Sfx{i};
  Sd((7:7)+(i-1)*7+19*NB, (7:7)+(i-1)*7+19*NB) = Sd2q{i};

  d((1:26)+(i-1)*26, 1) = [a( i,1:6), fB(i,1:6), f( i,1:6), tau(i,1), fx(i,1:6), d2q(i,1)]';
end



