function  S = scalc( jtyp )

% scalc motion subspace matrices.
% S=jcalc(type) returns the motion subspace matrix for a joint of the given
% type.  jtyp can be either a string or a structure containing a 
% string-valued field called 'code'.  Either way, the string contains the
% joint type code.  For joints that take parameters (e.g. helical), 
% jtyp must be a structure, and it must contain a field called 'pars', 
% which in turn is a structure containing one or more parameters.  (For a 
% helical joint, pars contains a parameter called 'pitch'.) 

if ischar( jtyp )
  code = jtyp;
else
  code = jtyp.code;
end

switch code
  case 'Rx'				% revolute X axis
    S = [1;0;0;0;0;0];
  case 'Ry'				% revolute Y axis
    S = [0;1;0;0;0;0];
  case {'R','Rz'}			% revolute Z axis
    S = [0;0;1;0;0;0];
  case 'Px'				% prismatic X axis
    S = [0;0;0;1;0;0];
  case 'Py'				% prismatic Y axis
    S = [0;0;0;0;1;0];
  case {'P','Pz'}			% prismatic Z axis
    S = [0;0;0;0;0;1];
  case 'H'				% helical (Z axis)
    S = [0;0;1;0;0;jtyp.pars.pitch];
  case 'r'				% planar revolute
    S = [1;0;0];
  case 'px'				% planar prismatic X axis
    S = [0;1;0];
  case 'py'				% planar prismatic Y axis
    S = [0;0;1];
  otherwise
    error( 'unrecognised joint code ''%s''', code );
end
