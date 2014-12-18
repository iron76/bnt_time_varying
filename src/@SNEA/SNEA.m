% SNEA solves inverse dynamics with redundant measurements.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation: 
%             mySNEA = SNEA(myModel, mySens).
% Computations are then performed with the method d = mySNEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef SNEA < deterministicIDsolver
   properties
      sparsified = 0;
      S
   end
   
   methods
      function b = SNEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b = b@deterministicIDsolver(m,y);                 
         b = initSparseIndices(b); 

      end % SNEA
      
      obj = initSparseIndices(obj)
      
      function disp(b)
         % Display SNEA object
         disp@deterministicIDsolver(b)
         fprintf('SNEA disp to be implemented! \n')
      end % disp
     
   end % methods
end % classdef

