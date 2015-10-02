function [muD, SD] = computePrior(obj)

NB = obj.IDmodel.modelParams.NB;
D = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
b = sparse(obj.ibs, ones(size(obj.ibs)), obj.bs, 19*NB, 1);

Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;

S_Dinv = Sv_inv;
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
bD     = b;

SD     = inv(D'*S_Dinv*D + S_dinv );
muD    = SD*( - D'*S_Dinv*bD );

muD  = muD(obj.id,1);
SD = SD(obj.id,obj.id);
