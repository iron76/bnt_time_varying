% LURNEA solves inverse dynamics with matrix form RNEA.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensRNEA(myModel)))
% with the following instanatiation: 
%             myLURNEA = LURNEA(myModel, mySens).
% Computations are then performed with the method d = myLURNEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef LUABA < deterministicMAPsolver
      
   methods
      function b = LUABA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b    = b@deterministicMAPsolver(m,y);
         
      end % LURNEA
     
   end % methods
end % classdef

