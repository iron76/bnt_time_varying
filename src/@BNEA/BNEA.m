% BNEA solves inverse dynamics with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB))
% and to describe the sensor distribution, e.g.:
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation:
%             myBNEA = SNEA(myModel, mySens).
% Computations are then performed with the method d = myBNEA.solveID()
% which relies on the Bayesian Network Toolbox (BNT) by Kevin Murphy.
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef BNEA < stochasticIDsolver
   
   properties
      bnt;
   end
   methods
      function b = BNEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b     = b@stochasticIDsolver(m,y);
         b.Sd  = zeros(b.IDmodel.modelParams.NB, b.IDmodel.modelParams.NB);
         b.bnt = b.bnetFromModel(b.IDmodel.modelParams, b.IDsens.sensorsParams);
      end % BNEA
      
      function disp(b)
         % Display BNEA object
         disp@deterministicIDsolver(b)
         fprintf('BNEA disp to be implemented! \n')
      end % disp
      
   end % methods
   
   methods(Static = true)
      bnt = bnetFromModel( dmodel , ymodel )
   end
   
end % classdef

