function [torquesM, baseReactionForceM ] = iDynTreeID(dyntree, qM, dqM, ddqM)
%IDYNTREEID Summary of this function goes here
%   Detailed explanation goes here

% set gravity
grav = iDynTree.SpatialAcc();
grav.setVal(2,-9.81);

dofs = dyntree.getNrOfDegreesOfFreedom();
q = iDynTree.VectorDynSize(dofs);
dq = iDynTree.VectorDynSize(dofs);
ddq = iDynTree.VectorDynSize(dofs);

torques = iDynTree.VectorDynSize(dofs);
baseReactionForce = iDynTree.Wrench();
    
torquesM =zeros(size(qM,1),dofs);
baseReactionForceM =zeros(size(qM,1),6);

for i = 1: size(qM,1)
    q.fromMatlab(qM(i,:));
    dq.fromMatlab(dqM(i,:));
    ddq.fromMatlab(ddqM(i,:));
    
    dyntree.setRobotState(q,dq,ddq,grav);
    
    % compute id with inverse dynamics
    dyntree.inverseDynamics(torques,baseReactionForce);
    
    torquesM(i,:) = torques.toMatlab();
    baseReactionForceM(i,:) = baseReactionForce.toMatlab();
end

end

