
function [I_O] = createSpatialInertia (I_cBar, m, c)


%CREATESPATIALINERTIA compute spatial inertia in Featherstone-like
%notatione (see Featherstone(2008),Rigid Body Dynamics Algorithms (2008), pg22).
%
%     I_0  spatial inertia of the rigid body with respect to the point O;

%     I_cBar   inertia of the rigid body with respect to the center of
%          mass C;
%     m    mass of the rigid body;
%     c    position of point 0 (the same one in URDF
%          <inertial> section).


I_O = [ I_cBar + m*skew(c)*(skew(c))'   m*skew(c);
              m*(skew(c))'              m*eye(3) ];



%TO DO :tests!



end



