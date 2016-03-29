function [] = priorAnalysis_plotResults(myBNEA,ll,bnetHat,i_learn)

    %clear 
    %close all
    %clc

    plotType = 'bar'; % or hinton;
    
    if(nargin<4)
        load('./experiments/humanFixedBase/intermediateDataFiles/savedBNet.mat');
        load('./experiments/humanFixedBase/intermediateDataFiles/EMResult.mat','myModel','mySens','bnet','bnetHat','ll');
    end

    if(strcmp(plotType, 'bar')==1)
            cov_ini = cell(length(i_learn),1);
            dispVect = i_learn;
            for j = 1:length(dispVect)
                cov_ini{j} = get_field(myBNEA.bnt.bnet.CPD{dispVect(j)}, 'cov');
                cov_ini_diag{j} = diag(get_field(myBNEA.bnt.bnet.CPD{dispVect(j)}, 'cov'));
            end

            fprintf('Initial y_fts covariance :');
            disp(cov_ini_diag{j}');
            fprintf('\n EM computed y_fts covariance \n');

            for j = 1:length(dispVect)
                s = size(cov_ini{j});
                cov_est_j = cov_ini{j};
                cov_est_j_diag = cov_ini_diag{j};

                for i = 1:length(ll)
                   loc_cov = get_field(bnetHat(i).CPD{dispVect(j)},'cov');
                   loc_cov_diag = diag(loc_cov);
                   cov_est_j = [cov_est_j loc_cov];
                   cov_est_j_diag = [cov_est_j_diag loc_cov_diag];

                   %% printing FT cov
                   if(dispVect(j)==17)
                       fprintf('Step %d, FT Cov diag : ',i);
                       disp(loc_cov_diag');
                   end
                end
                cov_est{j} = cov_est_j;
                cov_est_diag{j} = cov_est_j_diag;
           end

           for i = 1:length(dispVect)
                figure(i);
                fprintf('Node %d, cov dim (%d, %d)\n',i,size(cov_est{i}));
                bar3(1:1+length(ll),cov_est_diag{i}');
                title(myBNEA.bnt.nodes.labels{dispVect(i)});
           end
    else

        dir_ind = cell2mat(myBNEA.bnt.nodes.index);
        figure(1);
        plot(ll);
        %axis tight;
        xlabel('Epoch');
        ylabel('Log likelihood');
        grid on;


        for i = 1:18 %length(dir_ind(1:NB*6))+1 : length(dir_ind)
               cov_ini = get_field(myBNEA.bnt.bnet.CPD{dir_ind(i)}, 'cov');
               cov_est = get_field(bnetHat.CPD{dir_ind(i)},         'cov');
               cov_upd = cov_est - cov_ini;
               %fprintf('[INFO] %s was updated by %f \n', myBNEA.bnt.nodes.labels{dir_ind(i)}, norm(cov_upd)./norm(cov_ini));
               fprintf('[INFO] %s covariance of %f recomputed as %f \n', myBNEA.bnt.nodes.labels{dir_ind(i)}, norm(cov_ini),norm(cov_upd));

               %if(isempty(strfind(myMAP.IDsens.sensorsParams.labels{i},'d2q')))
               if(strcmp(diagType,'hinton') == 1) 
                   figure(1+ceil(i/2));
            %              if(numSubplots == 4)
            %                  subplot(2,2,mod(i,4)+1);
            %              elseif(numSubplots == 2)

                    subplot(2,2,mod(i+1,2)+1);
                    hintonDiagram2(cov_ini,'r');
                    title(strcat(myBNEA.bnt.nodes.labels{i},'-ini'));

                    subplot(2,2,mod(i+1,2)+3);
                    hintonDiagram2(cov_upd,'b');
                    title(strcat(myBNEA.bnt.nodes.labels{i},'-upd'));
               else
                   figure(ceil(i/4));
                     subplot(2,2,mod(i,4)+1);
                     bar([diag(cov_ini) diag(cov_est)]);
                     title( myBNEA.bnt.nodes.labels{i});
               end
                %bar([diag(cov_ini) diag(cov_est)]);

                     drawnow;

                % end
        end
    end
end