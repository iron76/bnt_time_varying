function [] = plotEvidence(t, evdSens, evdObt)

    str1 = 'sample';
    str2 = 'sampleTest';
    
    figure;
    subplot(2,2,1);
    plot(t,evdSens(21:23,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' mu']);
    subplot(2,2,2);
    plot(t,evdSens(24:26,:));axis tight;grid on;legend('x','y','z');
    title([str1,' F']);
    subplot(2,2,3);
    plot(t,evdObt(21:23,:));axis tight;grid on;legend('x','y','z');
    title([str2,' mu']);
    subplot(2,2,4);
    plot(t,evdObt(24:26,:));axis tight;grid on;legend('x','y','z');
    title([str2,' F']);
    
    %%%%%%%%%%%%
    figure;
    subplot(2,2,1);
    plot(t,evdSens(9:11,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' IMU-ang']);
    subplot(2,2,2);
    plot(t,evdObt(9:11,:));axis tight;grid on;legend('x','y','z');
    title([str2,' IMU-ang']);
    subplot(2,2,3);
    plot(t,evdSens(12:14,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' IMU-lin']);
    subplot(2,2,4);
    plot(t,evdObt(12:14,:));axis tight;grid on;legend('x','y','z');
    title([str2,' IMU-lin']);
    
    %%%%%%%%%%%%%
    figure;
    subplot(2,2,1);
    plot(t,evdSens([15:20,2:7],:)); axis tight;grid on;legend('x','y','z');
    title([str1,' Ftx']);
    subplot(2,2,2);
    plot(t,evdObt([15:20,2:7],:));axis tight;grid on;legend('x','y','z');
    title([str2,' Ftx']);
    subplot(2,2,3);
    plot(t,evdSens([1,8],:)); axis tight;grid on;legend('ddq1','ddq2');
    title([str1,' ddq']);
    subplot(2,2,4);
    plot(t,evdObt([1,8],:));axis tight;grid on;legend('ddq1','ddq2');
    title([str2,' ddq']);
    
end