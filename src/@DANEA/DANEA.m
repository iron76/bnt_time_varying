% DANEA solves inverse dynamics with accelerometers.
%
% The class is instantiated with two additional classes to describe the
% articulated rigid body model, e.g.:
%            myModel = model(autoTree(NB)) 
% and to describe the sensor distribution, e.g.: 
%            mySens  = sensors(autoSensStochastic(autoSensSNEA(myModel)))
% with the following instanatiation: 
%             myDANEA = DANEA(myModel, mySens).
% Computations are then performed with the method d = myDANEA.solveID();
%
% Author: Francesco Nori
% Genova, Dec 2014

classdef DANEA < derivativeACsolver
      
   methods
      function b = DANEA(m,y)
         if nargin == 0
            error(['You should provide a ' ...
               'model to instantiate SNEA'] )
         else
            if ~checkModel(m.modelParams)
               error(['You should provide a featherstone-like ' ...
                  'model to instantiate SNEA'] )
            end
         end
         b    = b@derivativeACsolver(m,y);
         b.Sd = zeros(b.IDmodel.modelParams.NB*26, b.IDmodel.modelParams.NB*26);
         
      end % DANEA
      
      function disp(b)
         % Display DANEA object
         disp@deterministicIDsolver(b)
         fprintf('DANEA disp to be implemented! \n')
      end % disp
      
   end % methods
end % classdef

