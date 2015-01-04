function obj = updateEvd(obj)

NB  = obj.IDmodel.modelParams.NB;

iy = 1;
for i = 1 : obj.IDsens.sensorsParams.ny
   obj.evd{obj.bnt.nodes.index{NB*6 + i}} = obj.IDmeas.y(iy:(iy+obj.IDsens.sensorsParams.sizes{i,1}-1));
   iy = iy + obj.IDsens.sensorsParams.sizes{i,1};
end