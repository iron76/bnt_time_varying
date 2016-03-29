function [] = plotMeasurements(t,y1,str1,y2,str2)

    figure;
    subplot(2,2,1);
    plot(t,y1(1:3,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' mu']);
    subplot(2,2,2);
    plot(t,y1(4:6,:));axis tight;grid on;legend('x','y','z');
    title([str1,' F']);
    subplot(2,2,3);
    plot(t,y2(1:3,:));axis tight;grid on;legend('x','y','z');
    title([str2,' mu']);
    subplot(2,2,4);
    plot(t,y2(4:6,:));axis tight;grid on;legend('x','y','z');
    title([str2,' F']);
    %%%%%%%%%%%%
    figure;
    subplot(2,2,1);
    plot(t,y1(7:9,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' IMU-ang']);
    subplot(2,2,2);
    plot(t,y2(7:9,:));axis tight;grid on;legend('x','y','z');
    title([str2,' mu-ang']);
    subplot(2,2,3);
    plot(t,y1(10:12,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' IMU-lin']);
    subplot(2,2,4);
    plot(t,y2(10:12,:));axis tight;grid on;legend('x','y','z');
    title([str2,' mu-lin']);
    %%%%%%%%%%%%%
    figure;
    subplot(2,2,1);
    plot(t,y1(13:24,:)); axis tight;grid on;legend('x','y','z');
    title([str1,' Ftx']);
    subplot(2,2,2);
    plot(t,y2(13:24,:));axis tight;grid on;legend('x','y','z');
    title([str2,' Ftx']);
    subplot(2,2,3);
    plot(t,y1(25:26,:)); axis tight;grid on;legend('ddq1','ddq2');
    title([str1,' ddq']);
    subplot(2,2,4);
    plot(t,y2(25:26,:));axis tight;grid on;legend('ddq1','ddq2');
    title([str2,' ddq']);
end