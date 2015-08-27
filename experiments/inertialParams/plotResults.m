figure
hold on
markerSize = 12;
linew = 3;
fontSize = 12;
plot(2:20,result{1,1}(2:end,1),'^','color','red','markerSize',markerSize,'linewidth',linew)
plot(2:20,result{1,2}(2:end,1),'+','color','blue','markerSize',markerSize,'linewidth',linew)
plot(2:20,result{1,3}(2:end,1),'o','color','green','markerSize',markerSize,'linewidth',linew)
xlabel('Iterations','FontSize',fontSize,'FontWeight','bold')
ylabel('Norm of errof of Base Parameters','FontSize',fontSize,'FontWeight','bold')
legend('JA','JA + GRF','JA + GRF + LA')
set(gca,'FontSize',9)
print('errorPlot','-depsc');

figure
plot(2:20,result{1,1}(2:end,2),'^','color','red','markerSize',markerSize,'linewidth',linew)
hold on
plot(2:20,result{1,2}(2:end,2),'+','color','blue','markerSize',markerSize,'linewidth',linew)
plot(2:20,result{1,3}(2:end,2),'o','color','green','markerSize',markerSize,'linewidth',linew)
xlabel('Iterations','FontSize',fontSize,'FontWeight','bold')
ylabel('Norm of errof of Not Identifiable Parameters','FontSize',fontSize,'FontWeight','bold')
legend('JA','JA + GRF','JA + GRF + LA')
set(gca,'FontSize',9)
print('errorNIPlot','-depsc');

