function e = residuals(obj,d)


pt(obj.id) = 1:length(obj.id);
d          = d(pt);

NB = obj.IDmodel.modelParams.NB;
D = sparse(obj.iDs, obj.jDs, obj.Ds, 19*NB, 26*NB);
b = sparse(obj.ibs, ones(size(obj.ibs)), obj.bs, 19*NB, 1);

% Dx = D(1:19*NB, 1:19*NB);
% Dy = D(1:19*NB, 19*NB+1:26*NB);

% Sv_inv = eye(19*NB)./sModel;
Sv_inv = obj.IDmodel.modelParams.Sv_inv.matrix;
% Sw_inv = eye(7*NB) ./sUknown;
Sw_inv = obj.IDmodel.modelParams.Sw_inv.matrix;
% Sw     = obj.IDmodel.modelParams.Sw.matrix;
% Sy_inv = eye(my)   ./sMeas;
Sy_inv = obj.IDsens.sensorsParams.Sy_inv.matrix;

Y = obj.IDsens.sensorsParams.Ys;

y      = obj.IDmeas.y;
S_Dinv = Sv_inv;
S_dinv = blkdiag(zeros(size(Sv_inv)), Sw_inv);
S_Yinv = Sy_inv;
bY     = zeros(size(y));
bD     = b;
muD    = zeros(length(S_dinv), 1);

e  = norm((D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y)*d-(Y'*S_Yinv*(y-bY) - D'*S_Dinv*bD + S_dinv * muD));
%norm(D'*S_Dinv*D + S_dinv + Y'*S_Yinv*Y - [D', Y']*[D;Y])
%norm(Y'*S_Yinv*(y-bY) - D'*S_Dinv*bD + S_dinv * muD - [D', Y']*[-bD; y])
return




