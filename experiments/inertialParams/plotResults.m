figure
hold on
markerSize = 8;
linew = 2;
fontSize = 12;

leng = size(result{1,2}(:,1),1);
plot(2:leng,result{1,2}(2:end,1),'+','color','blue','markerSize',markerSize,'linewidth',linew)
leng = size(result{1,3}(:,1),1);
plot(2:leng,result{1,3}(2:end,1),'o','color','green','markerSize',markerSize,'linewidth',linew)
leng = size(result{1,1}(:,1),1);
plot(2:leng,result{1,1}(2:end,1),'^','color','red','markerSize',markerSize,'linewidth',linew)

set(gca,'FontSize',fontSize,'FontWeight','bold')

xlabel('Iterations','FontSize',fontSize,'FontWeight','bold')
ylabel('Norm of errof of Base Parameters','FontSize',fontSize,'FontWeight','bold')
legend('JA + GRF','JA + GRF + LA','JA + GRF + LA + FT')
%set(gca,'FontSize',9)
print('errorPlot','-depsc');
print('errorPlot','-dpng');

figure


hold on
leng = size(result{1,2}(:,1),1);
plot(2:leng,result{1,2}(2:end,2),'+','color','blue','markerSize',markerSize,'linewidth',linew)

leng = size(result{1,3}(:,1),1);
plot(2:leng,result{1,3}(2:end,2),'o','color','green','markerSize',markerSize,'linewidth',linew)

leng = size(result{1,1}(:,1),1);
plot(2:leng,result{1,1}(2:end,2),'^','color','red','markerSize',markerSize,'linewidth',linew)


set(gca,'FontSize',fontSize,'FontWeight','bold')

xlabel('Iterations','FontSize',fontSize,'FontWeight','bold')
ylabel('Norm of errof of Not Identifiable Parameters','FontSize',fontSize,'FontWeight','bold')
legend('JA + GRF','JA + GRF + LA','JA + GRF + LA + FT')
%set(gca,'FontSize',9)
print('errorNIPlot','-depsc');
print('errorNIPlot','-dpng');

