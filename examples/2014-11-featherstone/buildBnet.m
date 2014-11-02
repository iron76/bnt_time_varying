function [ IDnet ] = buildBnet( model )

NB   = model.NB;
N    = NB*6;
sUkn = 1e12;

for i = 1 : NB
  nodes.index{i}    = (i-1)*6 + 1 : i*6;
  nodes.sizes{i,1}  = [6 6 6 1 6 1]';
  nodes.labels{i,1} = ['a' num2str(i)];
  nodes.labels{i,2} = ['fB' num2str(i)];
  nodes.labels{i,3} = ['f' num2str(i)];
  nodes.labels{i,4} = ['tau' num2str(i)];
  nodes.labels{i,5} = ['fx' num2str(i)];
  nodes.labels{i,6} = ['d2q' num2str(i)];
end

dag = zeros(N, N);
for i = 1 : NB
  % a{i} = ... + S{i}*qdd(i);
  dag(nodes.index{i}(6), nodes.index{i}(1)) = 1;
  % fB{i} = model.I{i}*a{i} + ... ;
  dag(nodes.index{i}(1), nodes.index{i}(2)) = 1;
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  dag(nodes.index{i}(2), nodes.index{i}(3)) = 1;
  dag(nodes.index{i}(5), nodes.index{i}(3)) = 1;
  % tau(i,1) = S{i}' * f{i};
  dag(nodes.index{i}(3), nodes.index{i}(4)) = 1;
  
  
  if model.parent(i) ~= 0
    % a{i} = Xup{i}*a{model.parent(i)} + ...
    dag(nodes.index{model.parent(i)}(1), nodes.index{i}(1)) = 1;
  end
  
  ind_j = find(model.parent == i);
  for j = ind_j
    % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
    dag(nodes.index{j}(3), nodes.index{i}(3)) = 1;
  end
end

dnodes = [];  % no discrete nodes
ns     = cell2mat(nodes.sizes);
bnet   = mk_bnet(dag, ns, 'discrete', dnodes);

link = zeros(N,1);
for i = 1 : NB
  for j = 1 : 6
    nodes.index{i}(j) = find(bnet.order == nodes.index{i}(j));
    link(nodes.index{i}(j), 1) = i;
  end
end
 

identity = eye(N);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);
bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);

for i = 1 : NB
  fxi  = nodes.index{i}(5);
  nfxi = nodes.sizes{i,1}(5);
  bnet.CPD{fxi} = gaussian_CPD(bnet, fxi, 'mean', zeros(nfxi,1), 'cov', sUkn.*eye(nfxi));
  
  d2qi  = nodes.index{i}(6);
  nd2qi = nodes.sizes{i,1}(6);
  bnet.CPD{d2qi} = gaussian_CPD(bnet, d2qi, 'mean', zeros(nd2qi,1), 'cov', sUkn.*eye(nd2qi));
end

IDnet.bnet    = bnet;
IDnet.nodes   = nodes;
IDnet.link    = link;
IDnet.engine  = [];
end

