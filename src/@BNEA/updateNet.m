function obj = updateNet(obj)

NB     = obj.IDmodel.modelParams.NB;

for i = 1:obj.IDmodel.modelParams.NB
   if obj.IDmodel.modelParams.parent(i) == 0
      % a{i} = Xup{i}*(-a_grav) + S{i}*qdd(i);
      ai   = obj.bnt.nodes.index{i}(1);
      nai  = obj.bnt.nodes.sizes{i,1}(1);
      d2qi = obj.bnt.nodes.index{i}(6);
      pars  = cell2mat(obj.bnt.bnet.parents(ai));
      Wa    = cell(1, length(pars));
      for j = 1 : length(pars)
         if pars(j) == d2qi
            Wa{1, j} = obj.IDmodel.S{i};
         end
      end
      W = cell2mat(Wa);
      obj.bnt.bnet.CPD{ai} = gaussian_CPD(obj.bnt.bnet, ai, 'mean', obj.Xup{i}*(-obj.IDmodel.g), 'cov', obj.IDmodel.modelParams.Sm.a{i}, 'weights', W);
   else
      % a{i} = ... + S{i}*qdd(i) + crm(v{i})*vJ;
      % a{i} = Xup{i}*a{model.parent(i)} + ...
      ai      = obj.bnt.nodes.index{i}(1);
      d2qi    = obj.bnt.nodes.index{i}(6);
      nai     = obj.bnt.nodes.sizes{i,1}(1);
      aj      = obj.bnt.nodes.index{obj.IDmodel.modelParams.parent(i)}(1);
      
      pars  = cell2mat(obj.bnt.bnet.parents(ai));
      Wa    = cell(1, length(pars));
      for j = 1 : length(pars)
         if pars(j) == aj
            Wa{1, j} = obj.Xup{i};
         elseif pars(j) == d2qi
            Wa{1, j} = obj.IDmodel.S{i};
         end
      end
      W = cell2mat(Wa);
      obj.bnt.bnet.CPD{ai} = gaussian_CPD(obj.bnt.bnet, ai, 'mean', crm(obj.v(:,i))*obj.vJ(:,i), 'cov', obj.IDmodel.modelParams.Sm.a{i}, 'weights', W);
   end
   % fB{i} = model.I{i}*a{i} + crf(v{i})*model.I{i}*v{i};
   fBi  = obj.bnt.nodes.index{i}(2);
   nfBi = obj.bnt.nodes.sizes{i,1}(2);
   pars  = cell2mat(obj.bnt.bnet.parents(fBi));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == ai
         Wa{1, j} = obj.IDmodel.modelParams.I{i};
      end
   end
   W = cell2mat(Wa);
   obj.bnt.bnet.CPD{fBi} = gaussian_CPD(obj.bnt.bnet, fBi, 'mean', crf(obj.v(:,i))*obj.IDmodel.modelParams.I{i}*obj.v(:,i), 'cov', obj.IDmodel.modelParams.Sm.fB{i}, 'weights', W);
end

for i = obj.IDmodel.modelParams.NB:-1:1
   % f{i} = fB{i} - Xa{i}' \ f_ext{i};
   % f{model.parent(j)} = f{model.parent(j)} + Xup{j}'*f{j};
   fBi = obj.bnt.nodes.index{i}(2);
   fi  = obj.bnt.nodes.index{i}(3);
   fxi = obj.bnt.nodes.index{i}(5);
   nfi = obj.bnt.nodes.sizes{i,1}(3);
   pars  = cell2mat(obj.bnt.bnet.parents(fi));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == fBi
         Wa{1, j} = eye(nfi);
      elseif pars(j) == fxi
         Wa{1, j} = -inv(obj.Xa{i}');
      else
         Wa{1, j} = obj.Xup{obj.bnt.link(pars(j))}';
      end
   end
   W = cell2mat(Wa);
   
   obj.bnt.bnet.CPD{fi} = gaussian_CPD(obj.bnt.bnet, fi, 'mean', zeros(nfi,1), 'cov', obj.IDmodel.modelParams.Sm.f{i}, 'weights', W);
   
   % tau(i,1) = S{i}' * f{i};
   taui  = obj.bnt.nodes.index{i}(4);
   ntaui = obj.bnt.nodes.sizes{i,1}(4);
   pars  = cell2mat(obj.bnt.bnet.parents(taui));
   Wa    = cell(1, length(pars));
   for j = 1 : length(pars)
      if pars(j) == fi
         Wa{1, j} = obj.IDmodel.S{i}';
      end
   end
   W = cell2mat(Wa);
   obj.bnt.bnet.CPD{taui} = gaussian_CPD(obj.bnt.bnet, taui, 'mean', zeros(ntaui, 1), 'cov', obj.IDmodel.modelParams.Sm.tau{i}, 'weights', W);
   
   % fxi
   fxi  = obj.bnt.nodes.index{i}(5);
   nfxi = obj.bnt.nodes.sizes{i,1}(5);
   obj.bnt.bnet.CPD{fxi} = gaussian_CPD(obj.bnt.bnet, fxi, 'mean', zeros(nfxi, 1), 'cov', obj.IDmodel.modelParams.Su.fx{i});
   
   % d2qi
   d2qi  = obj.bnt.nodes.index{i}(6);
   nd2qi = obj.bnt.nodes.sizes{i,1}(6);
   obj.bnt.bnet.CPD{d2qi} = gaussian_CPD(obj.bnt.bnet, d2qi, 'mean', zeros(nd2qi, 1), 'cov', obj.IDmodel.modelParams.Su.d2q{i});
   
end

for i = 1 : obj.IDsens.sensorsParams.ny
   yi  = obj.bnt.nodes.index{NB*6 + i};
   nyi = obj.bnt.nodes.sizes{NB*6 + i};
   obj.bnt.bnet.CPD{yi} = gaussian_CPD(obj.bnt.bnet, yi, 'mean', zeros(nyi, 1), 'cov', obj.IDsens.sensorsParams.Sy{i,1}, 'weights', obj.bnt.Wy{i,1});
end