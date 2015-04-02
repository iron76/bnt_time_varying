function  [Xj,S,dXjdq] = jcalcderiv( jtyp, q )

% jcalcderiv computes joint transform, motion subspace matrices and joint
% transform derivatives. Xj and S are the same output of jcalc
% documentation, while dXjdq is the derivarive of Xj with respect to the
% state position variable q. 
% [Xj,S]=jcalcderiv(type,q)  
if ischar( jtyp )
  code = jtyp;
else
  code = jtyp.code;
end

switch code
  case 'Rx'				% revolute X axis
    Xj = rotx(q);
    S = [1;0;0;0;0;0];
    dXjdq = rotxderiv(q);
  case 'Ry'				% revolute Y axis
    Xj = roty(q);
    S = [0;1;0;0;0;0];
    dXjdq = rotyderiv(q);
  case {'R','Rz'}			% revolute Z axis
    Xj = rotz(q);
    S = [0;0;1;0;0;0];
    dXjdq = rotzderiv(q);
  case 'Px'				% prismatic X axis
    Xj = xlt([q 0 0]);
    S = [0;0;0;1;0;0];
    dXjdq = traslxderiv(q);
  case 'Py'				% prismatic Y axis
    Xj = xlt([0 q 0]);
    S = [0;0;0;0;1;0];
    dXjdq = traslyderiv(q);
  case {'P','Pz'}			% prismatic Z axis
    Xj = xlt([0 0 q]);
    S = [0;0;0;0;0;1];
    dXjdq = traslzderiv(q);
  case 'H'				% helical (Z axis)
    Xj = rotz(q) * xlt([0 0 q*jtyp.pars.pitch]);
    S = [0;0;1;0;0;jtyp.pars.pitch];
    error( 'joint code ''%s'' not supported by jcalcderiv', code );
  case 'r'				% planar revolute
    Xj = plnr( q, [0 0] );
    S = [1;0;0];
    error( 'joint code ''%s'' not supported by jcalcderiv', code );
  case 'px'				% planar prismatic X axis
    Xj = plnr( 0, [q 0] );
    S = [0;1;0];
    error( 'joint code ''%s'' not supported by jcalcderiv', code );
  case 'py'				% planar prismatic Y axis
    Xj = plnr( 0, [0 q] );
    S = [0;0;1];
    error( 'joint code ''%s'' not supported by jcalcderiv', code );
  otherwise
    error( 'unrecognised joint code ''%s''', code );
end
