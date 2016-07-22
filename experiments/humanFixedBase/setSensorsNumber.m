function [ myMAP, Ymatrix,b_Y,yMeas ] = setSensorsNumber(dmodel, myModel,len,q,dq,ddq,currentParams,currentTrialSens, data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



        %% Build sensor model                  
        
        sens.parts    = {'torso'};      
        sens.labels   = {'imu'};
        sens.ndof     = {6};
% %         sens.parts    = {'leg'         ,'torso'};       %force of the forceplate is trasmitted into the leg
% %         sens.labels   = {'fts'         ,'imu'  };
% %         sens.ndof     = {6             ,6      };
% %          label_to_plot = {'fts'         ,'imu'  };
        ymodel  = humanThreeLinkSens(dmodel, sens);  
        ymodel  = humanThreeLinkSensStochastic(ymodel);     
        mySens  = sensors(ymodel);  
        myMAP   = MAP(myModel, mySens);

         %% Build data.y anda data.Sy 
        % data.y are ordered in:  
        % - angular-linear notation 
        % - the form [a2 ftx1 ftx2 ddq1 ddq2]
        % - sensor frame

        %===== data.y
        data.y  = [];
        for i = 1 : length(sens.labels)
             eval(['data.y  = [data.y  data.ys_sensFrame_' sens.labels{i} '];']);
        end

        % Add null external forces ftx = 0
        data.y  = [data.y, zeros(len,6*dmodel.NB)];
                
        % Add ddq measurements
        data.y  = [data.y, data.ddq];

        yMeas = data.y';
        
        %% Build Ymatrix manually
        % Y has to be consistent with measurements form [ddq1 ddq2]

        Y = zeros(myMAP.IDsens.m, 26*dmodel.NB);
        Y(1:6,20:25) = eye(6);
        Y(7,26) = eye(1);
        Y(8:13,27:32) = currentTrialSens.X_imu_2;
        Y(14:19,46:51) = eye(6);
        Y(20,52) = eye(1);
        
% %         Y_4sens = cell2mat(ymodel_4sens.Y);
% %         %Y_4sens = zeros (ymodel.m,26*dmodel.NB);
% %         Y_4sens(7:12,27:32) % %         Y_4sens = cell2mat(ymodel_4sens.Y);
% %         %Y_4sens = zeros (ymodel.m,26*dmodel.NB);
% %         Y_4sens(7:12,27:32) = X_imu_2;

% %         Y_4sens(13:18,20:25) = eye(6);
% %         Y_4sens(19:24,46:51) = eye(6);
% %         Y_4sens(25,26) = eye(1);
% %         Y_4sens(26,52) = eye(1);


        Ymatrix = cell(len,1);
        for i = 1 : len
            %the only row in Ymatrix that is time varying
            %Y(1:6,13:18) = XStar_fp_0 * XStar_0_1{i};
            Ymatrix{i} = Y;
        end
        
         %% Computing v
         
        ymodel_RNEA  = autoSensRNEA(dmodel);
        mySens_RNEA  = sensors(ymodel_RNEA);
        myRNEA       = RNEA(myModel, mySens_RNEA);

        y_RNEA_f = zeros(6*dmodel.NB, len);
        y_RNEA_ddq = zeros(dmodel.NB, len);
        fx = cell(dmodel.NB);

        %Ordering y_RNEA in the form [fx1 fx2 ddq1 ddq2]
        for i = 1 : dmodel.NB
            for j = 1 : len
                fx{i,1} = zeros(6,1); 
                y_RNEA_f(((1:6)+(i-1)*6), j) = [fx{i,1}];
                y_RNEA_ddq(i, j) = [ddq(i,j)];
            end
            y_RNEA1 = [y_RNEA_f ; y_RNEA_ddq];
            y_RNEA2 = [y_RNEA_f ; y_RNEA_ddq];
        end

        d = zeros (26*myRNEA.IDmodel.modelParams.NB,len);
        v_RNEA = cell(len,1);
        resRNEA.tau = zeros(length(len),dmodel.NB);
        
        for i = 1 : len
             myRNEA = myRNEA.setState(q(:,i), dq(:,i));
             myRNEA = myRNEA.setY(y_RNEA2(:,i));
             myRNEA = myRNEA.solveID();

             d(:,i) = myRNEA.d;
             
             v_RNEA{i,:} = myRNEA.v; 
             resRNEA.tau(i,:) = myRNEA.tau;
        end
        
        clear fx;
        
        %% Build bias b_Y manually
        % b_Y has to be consistent with Ymatrix

%         footMass =  currentParams.footMass;
%         posP_0 = [0; 0; (0.5*currentParams.footHeight)];
%         footIxx =  currentParams.footIxx;
%         footIyy =  currentParams.footIyy;
%         footIzz =  currentParams.footIzz;

        b_Y = zeros (size(data.y')); 
        R_imu_2 = currentTrialSens.X_imu_2(1:3,1:3);

%         a_G_grav = [0;0;0;0;0;-9.8100]; %Featherstone-like notation, in global reference
%         X_G_0 = currentTrialSens.X_G_0;
%         X_0_G = InverseAdjTransform(X_G_0);
%         a_0_grav = X_0_G * a_G_grav;

%         I_0 = createSpatialInertia(footIxx,footIyy,footIzz,footMass,posP_0);

        %b_Y(1:6,1:len)   = repmat((-XStar_fp_0 * I_0 * a_0_grav),1,len);
        
%         v = cell(size(q))';
%         fx = zeros (6,1);
%         fext    = cell(1,2);
%         for i = 1 : dmodel.NB
%              fext{i}    = fx;
%         end

        % exploiting velocitiy v_RNEA coming from RNEA class computation
        for i = 1 : len 
            A =R_imu_2*v_RNEA{i,1}(1:3,2);
            B =((currentTrialSens.X_imu_2(4:6,1:3)*v_RNEA{i,1}(1:3,2))+(R_imu_2*v_RNEA{i,1}(4:6,2)));
            b_Y(4:6,i) = cross(A,B);
        end 
        
        clear A;
        clear B;
%         clear footIxx;
%         clear footIyy;
%         clear footIzz;
%         
        
end

