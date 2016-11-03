% clear variables and close windows
clear all
close all
clc

step=10;

x=(1:8*step)/step;
y1=gampdf(x,1,2);
y2=gampdf(x,2,1);
y3=gampdf(x,3,0.5);

figure(1)
plot(x,y1,'k','LineWidth',2);
hold on
plot(x,y2,':r','LineWidth',2);
plot(x,y3,'--','LineWidth',2);
xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
set(gca,'Ytick',[0.0:0.1:0.6],'YTickLabel',[0.0:0.1:0.6],'FontSize',16,'FontWeight','Bold')
set(gca,'Xtick',[0:2:8],'XTickLabel',{0,2,4,6,8},'FontSize',16,'FontWeight','Bold')
title('Gamma densities','FontSize',16,'FontWeight','Bold');
ylim([-0.01 0.6]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on
% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss04_01.pdf
%print -painters -dpng -r600 STFloss04_01.png

figure(2)
semilogy(x,y1,'k','LineWidth',2);
hold on
semilogy(x,y2,':r','LineWidth',2);
semilogy(x,y3,'--','LineWidth',2);

ylim([10e-4, 10e-1]);
set(gca,'Xtick',[0:2:8],'XTickLabel',{0,2,4,6,8},'FontSize',16,'FontWeight','Bold')
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Gamma densities','FontSize',16,'FontWeight','Bold');
% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss04_02.pdf
%print -painters -dpng -r600 STFloss04_02.png

figure(3)

x=(1:25*step)/step;
y1=lognpdf(x,2,1);
y2=lognpdf(x,2,0.1);
y3=lognpdf(x,0.5,2);

plot(x,y1,'k' ,'LineWidth',2);
hold on
plot(x,y2,':r' ,'LineWidth',2);
plot(x,y3,'--' ,'LineWidth',2);
xlim([0 25]);
ylim([-0.01 0.6])

xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Log-normal densities','FontSize',16,'FontWeight','Bold');
set(gca,'Ytick',[0.0:0.1:0.6],'YTickLabel',[0.0:0.1:0.6],'FontSize',16,'FontWeight','Bold')
set(gca,'Xtick',[0:5:25],'XTickLabel',[0:5:25],'FontSize',16,'FontWeight','Bold')
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on
% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss04_03.pdf
%print -painters -dpng -r600 STFloss04_03.png

figure(4)

semilogy(x,y1,'k' ,'LineWidth',2);
hold on
semilogy(x,y2,':r' ,'LineWidth',2);
semilogy(x,y3,'--' ,'LineWidth',2);
xlim([0 25]);
ylim([10e-4 10e-1])

xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Log-normal densities','FontSize',16,'FontWeight','Bold');
set(gca,'Xtick',[0:5:25],'XTickLabel',[0:5:25],'FontSize',16,'FontWeight','Bold')
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss04_04.pdf
%print -painters -dpng -r600 STFloss04_04.png