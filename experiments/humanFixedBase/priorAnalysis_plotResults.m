function [] = priorAnalysis_plotResults(myBNEA,ll,bnetHat,i_learn,mySens)

    close all
    if(nargin<5)
        %% loading from memory if not called as a function
        close all; 
        clear;
        fprintf('Loading from saved result\n');
        load('./experiments/humanFixedBase/intermediateDataFiles/savedBNet.mat');
        load('./experiments/humanFixedBase/intermediateDataFiles/EMResult.mat');
    end

    
    plotType = 'bar'; % or hinton;
    
    cov_ini = cell(length(i_learn),1);
    dispVect = i_learn;
    for j = 1:length(dispVect)
        cov_ini{j} = struct(myBNEA.bnt.bnet.CPD{dispVect(j)}).cov;
        cov_ini_diag{j} = diag(struct(myBNEA.bnt.bnet.CPD{dispVect(j)}).cov);
    end

    %fprintf('Initial y_fts covariance :');
    %disp(cov_ini_diag{j}');
   % fprintf('\n EM computed y_fts covariance \n');
    
    for j = 1:length(dispVect)
        s = size(cov_ini{j});
        cov_est_j = cov_ini{j};
        
        covStore{j}=zeros(length(ll)+1,size(cov_ini{j},1),size(cov_ini{j},2));
        covStore{j}(1,:,:) = cov_ini{j};
        
        cov_est_j_diag = cov_ini_diag{j};
        
        if(j == 3 || j == 4) %% IMU or FT
            temp = mySens.sensorsParams.sensXLink{j-2}*cov_ini{j}*mySens.sensorsParams.sensXLink{j-2}';
            cov_est_j_diag_sensFrame = diag(temp);
            
            if(j == 4)
                fprintf('Ini, y_FT Cov :\n');
                disp(temp);
            end
        end
        
        for i = 1:length(ll)
           loc_cov = struct(bnetHat(i).CPD{dispVect(j)}).cov;
           covStore{j}(i+1,:,:)=loc_cov;
           loc_cov_diag = diag(loc_cov);
           cov_est_j = [cov_est_j loc_cov];
           cov_est_j_diag = [cov_est_j_diag loc_cov_diag];
           
           if(j == 3 || j == 4)
            loc_cov_sensFrame =  mySens.sensorsParams.sensXLink{j-2}*loc_cov*mySens.sensorsParams.sensXLink{j-2}';
            cov_est_j_diag_sensFrame =  [cov_est_j_diag_sensFrame diag(loc_cov_sensFrame)];
            
            if(j == 4);
               fprintf('Run %d, y_FT Cov :\n',i);
               disp(loc_cov_sensFrame);               
            end
           end
           
           
        end
        cov_est{j} = cov_est_j;
        cov_est_diag{j} = cov_est_j_diag;
        
        if(j == 3 || j == 4) %% IMU or FT
           cov_est_diag_sensFrame{j} = cov_est_j_diag_sensFrame;
        end
    end

   %% bar3d plot of the main diagonals of the obtained covariances
   for i = 1:length(dispVect)
        lastFig = figure(i);
        if( i == 3 || i == 4)
            subplot(1,2,1);
        end
        fprintf('Node %d, cov dim (%d, %d)\n',i,size(cov_est{i}));
        bar3(1:1+length(ll),cov_est_diag{i}'); axis tight;
        title(sprintf('%s (linkFrame)',myBNEA.bnt.nodes.labels{dispVect(i)}));
        
        if( i == 3 || i == 4)
            subplot(1,2,2);
            bar3(1:1+length(ll),cov_est_diag_sensFrame{i}'); axis tight;
            title(sprintf('%s (sensorFrame)',myBNEA.bnt.nodes.labels{dispVect(i)}));
        end
        %title(myBNEA.bnt.nodes.labels{dispVect(i)});
   end

   %% covariance visulaization for IMU and FTS
   covStoreSenseFrame = covStore;
   for i = [3,4]%1:length(dispVect)
        for j = 1:4%length(ll)+1;
            figure(i+lastFig);
            subplot(2,2,j);
            imagesc(squeeze(covStore{i}(j,:,:)));
            set(gca,'xTickLabel',1:6);
            set(gca,'yTickLabel',1:6);

            title(sprintf('%s (linkFrame), run %d',myBNEA.bnt.nodes.labels{dispVect(i)},j-1));

          %  if(j == 3 || j == 4); %% IMU
             covStoreSenseFrame{i}(j,:,:) = mySens.sensorsParams.sensXLink{i-2}*squeeze(covStore{i}(j,:,:))*mySens.sensorsParams.sensXLink{i-2}';
           % end            
            
            figure(i+lastFig+2);
            subplot(2,2,j);
            imagesc(squeeze(covStoreSenseFrame{i}(j,:,:)));
            set(gca,'xTickLabel',1:6);
            set(gca,'yTickLabel',1:6);

            title(sprintf('%s (sensorFrame), run %d',myBNEA.bnt.nodes.labels{dispVect(i)},j-1));

        end
   end

   
figure;

%plot(-log(-ll),'-o','lineWidth',2.0); axis tight; hold on;
semilogy(ll,'-o','lineWidth',2.0); axis tight; hold on;
xlabel('EM run');
ylabel('Log(Log likelihood)');
title('EM propogation');
   
end
