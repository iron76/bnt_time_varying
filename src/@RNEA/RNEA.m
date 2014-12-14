% RNEA computes inverse dynamics with recursive Newton-Euler
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensRNEA(myModel)))
% with the following instanatiation: 
%             myRNEA = RNEA(myModel, mySens).
% Computations are then performed with the method d = myRNEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef RNEA < deterministicIDsolver
  
   methods
      function b = RNEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate RNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate RNEA'] )
            end
         end
         b = b@deterministicIDsolver(m,y);
      end % RNEA
      
      function disp(b)
         % Display RNEA object
         disp@deterministicIDsolver(b)
         fprintf('RNEA disp to be implemented! \n')
      end % disp
      

   end % methods
   
end % classdef

