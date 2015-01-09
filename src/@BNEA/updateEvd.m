function obj = updateEvd(obj)

NB  = obj.IDmodel.modelParams.NB;

iy = 1;
for i = 1 : obj.IDsens.sensorsParams.ny
   obj.evd{obj.bnt.nodes.index{NB + i}} = obj.IDmeas.y(iy:(iy+obj.IDsens.sensorsParams.sizes{i,1}-1));
   iy = iy + obj.IDsens.sensorsParams.sizes{i,1};
end

if strcmp(obj.inf_engine, 'jtree_inf_engine')
   %% jtree
   if isempty(obj.bnt.engine)
      obj.bnt.engine = jtree_inf_engine(obj.bnt.bnet);
      obj.bnt.engine = enter_evidence(obj.bnt.engine, obj.evd);
   else
      % disp('Using previous path on jtree!')
      obj.bnt.engine = bnet_to_engine(obj.bnt.bnet, obj.bnt.engine);
      obj.bnt.engine = enter_evidence(obj.bnt.engine, obj.evd);
   end
elseif strcmp(obj.inf_engine, 'gaussian_inf_engine')
   %%gaussian
   obj.bnt.engine = gaussian_inf_engine(obj.bnt.bnet);
   obj.bnt.engine = enter_evidence(obj.bnt.engine, obj.evd);
end