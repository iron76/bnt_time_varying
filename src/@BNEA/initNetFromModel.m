function [ IDnet ] = initNetFromModel( dmodel , ymodel )

NB        = dmodel.NB;
N         = NB*6+ymodel.ny;

for i = 1 : NB
  nodes.index{i}    = (i-1)*6 + 1 : i*6;
  nodes.sizes{i,1}  = [6 6 6 1 6 1]';
  nodes.labels{nodes.index{i}(1)} = ['a' num2str(i)];
  nodes.labels{nodes.index{i}(2)} = ['fB' num2str(i)];
  nodes.labels{nodes.index{i}(3)} = ['f' num2str(i)];
  nodes.labels{nodes.index{i}(4)} = ['tau' num2str(i)];
  nodes.labels{nodes.index{i}(5)} = ['fx' num2str(i)];
  nodes.labels{nodes.index{i}(6)} = ['d2q' num2str(i)];
end

for i = 1 : ymodel.ny
  nodes.index{i+NB}    = i+NB*6;
  nodes.sizes{i+NB,1}  = ymodel.sizes{i,1};
  nodes.labels{nodes.index{i+NB}} = ymodel.labels{i,1};
end

dag = zeros(N, N);
for i = 1 : NB
  % a{i} = ... + S{i}*qdd(i);
  dag(nodes.index{i}(6), nodes.index{i}(1)) = 1;
  % fB{i} = dmodel.I{i}*a{i} + ... ;
  dag(nodes.index{i}(1), nodes.index{i}(2)) = 1;
  % f{i} = fB{i} - Xa{i}' \ f_ext{i};
  dag(nodes.index{i}(2), nodes.index{i}(3)) = 1;
  dag(nodes.index{i}(5), nodes.index{i}(3)) = 1;
  % tau(i,1) = S{i}' * f{i};
  dag(nodes.index{i}(3), nodes.index{i}(4)) = 1;
  
  
  if dmodel.parent(i) ~= 0
    % a{i} = Xup{i}*a{dmodel.parent(i)} + ...
    dag(nodes.index{dmodel.parent(i)}(1), nodes.index{i}(1)) = 1;
  end
  
  ind_j = find(dmodel.parent == i);
  for j = ind_j
    % f{dmodel.parent(j)} = f{dmodel.parent(j)} + Xup{j}'*f{j};
    dag(nodes.index{j}(3), nodes.index{i}(3)) = 1;
  end
  
  for j = 1 : ymodel.ny
     if ~isempty(find(ymodel.Y{j,i}(:,1:6),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,1:6);
        dag(nodes.index{i}(1), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(1);
     end
     if ~isempty(find(ymodel.Y{j,i}(:,7:12),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,7:12);
        dag(nodes.index{i}(2), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(2);
     end
     if ~isempty(find(ymodel.Y{j,i}(:,13:18),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,13:18);
        dag(nodes.index{i}(3), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(3);
     end
     if ~isempty(find(ymodel.Y{j,i}(:,19:19),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,19:19);
        dag(nodes.index{i}(4), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(4);
     end
     if ~isempty(find(ymodel.Y{j,i}(:,20:25),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,20:25);
        dag(nodes.index{i}(5), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(5);
     end
     if ~isempty(find(ymodel.Y{j,i}(:,26:26),1))
        IDnet.Wy{j,1}     = ymodel.Y{j,i}(:,26:26);
        dag(nodes.index{i}(6), nodes.index{j+NB}) = 1;
        IDnet.child{j,1} = nodes.index{i}(6);
     end
  end
end

dnodes = [];  % no discrete nodes
ns     = cell2mat(nodes.sizes);
bnet   = mk_bnet(dag, ns, 'discrete', dnodes);

link = zeros(N,1);
for i = 1 : NB
  for h = 1 : 6
    nodes.index{i}(h) = find(bnet.order == nodes.index{i}(h));
    link(nodes.index{i}(h), 1) = i;
  end
end

for i = 1 : ymodel.ny
  nodes.index{NB + i} = find(bnet.order == nodes.index{NB + i});
end

identity = eye(N);
new_dag = dag*identity(:,bnet.order);
new_dag = identity(bnet.order, :)*new_dag;
new_ns = ns(bnet.order);
bnet = mk_bnet(new_dag, new_ns, 'discrete', dnodes);


IDnet.bnet   = bnet;
IDnet.nodes  = nodes;
IDnet.link   = link;
IDnet.engine = [];
end

