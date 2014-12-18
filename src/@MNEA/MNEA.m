% MNEA solves inverse dynamics with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation: 
%             myMNEA = SNEA(myModel, mySens).
% Computations are then performed with the method d = myMNEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef MNEA < stochasticIDsolver
      
   methods      
      function b = MNEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate MNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate MNEA'] )
            end
         end
         b    = b@stochasticIDsolver(m,y);
         b.Sd = zeros(b.IDmodel.modelParams.NB, b.IDmodel.modelParams.NB);         
      end % MNEA
      
      function disp(b)
         % Display SNEA object
         disp@deterministicIDsolver(b)
         fprintf('MNEA disp to be implemented! \n')
      end % disp
      
   end % methods
end % classdef

