% ANEA solves accelerations with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation: 
%             myANEA = ANEA(myModel, mySens).
% Computations are then performed with the method a = myANEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef ANEA < stochasticACsolver
      
   methods
      function b = ANEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b    = b@stochasticACsolver(m,y);
         b.Sd = zeros(b.IDmodel.modelParams.NB*7, b.IDmodel.modelParams.NB*7);
         
      end % ANEA
      
      function disp(b)
         % Display SNEA object
         disp@deterministicIDsolver(b)
         fprintf('ANEA disp to be implemented! \n')
      end % disp
      
   end % methods
end % classdef

