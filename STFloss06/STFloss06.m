% ---------------------------------------------------------------------
% clear variables and close windows
clear all
close all
clc


step=20;

% Burr densities
x=(1:144*step)/step;
y1=Burrpdf(x,0.5,2,1.5);
y2=Burrpdf(x,0.5,0.5,5);
y3=Burrpdf(x,2,1,0.5);

figure(1)
plot(x,y1,'k','LineWidth',2);
hold on
plot(x,y2,':r','LineWidth',2);
plot(x,y3,'--','LineWidth',2);
xlim([0,8])
xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Burr densities','FontSize',16,'FontWeight','Bold');
ylim([0 1.2]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on
% to save the plot in pdf or png please uncomment next 2 lines:
print -painters -dpdf -r600 STFloss06_01.pdf
print -painters -dpng -r600 STFloss06_01.png

% Weibull densities
x=(1:144*step)/step;
y1=wblpdf(x,1.^(-1./0.5),0.5);
y2=wblpdf(x,1.^(-1/2),2);
y3=wblpdf(x,0.01.^(-1/6),6);


figure(2)
plot(x,y1,'k','LineWidth',2);
hold on
plot(x,y2,':r','LineWidth',2);
plot(x,y3,'--','LineWidth',2);
xlim([0,5])
xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Weibull densities','FontSize',16,'FontWeight','Bold');
ylim([-0.01 1.2]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on
% to save the plot in pdf or png please uncomment next 2 lines:
print -painters -dpdf -r600 STFloss06_02.pdf
print -painters -dpng -r600 STFloss06_02.png

% Burr double-logarithmic densities
figure(3)

x=(1:144*step)/step;
y1=Burrpdf(x,0.5,2,1.5);
y2=Burrpdf(x,0.5,0.5,5);
y3=Burrpdf(x,2,1,0.5);

loglog(x,y1,'k','LineWidth',2);
hold on
loglog(x,y2,':r','LineWidth',2);
loglog(x,y3,'--','LineWidth',2);

xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Burr densities','FontSize',16,'FontWeight','Bold');
ylim([10e-5, 10e-1]);
xlim([10e-2, 10e1]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');

box on
% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss06_03.pdf
%print -painters -dpng -r600 STFloss06_03.png

% Weibull double-logarithmic densities

x=(1:144*step)/step;
y1=wblpdf(x,1.^(-1./0.5),0.5);
y2=wblpdf(x,1.^(-1/2),2);
y3=wblpdf(x,0.01.^(-1/6),6);


figure(4)
semilogy(x,y1,'k','LineWidth',2);
hold on
semilogy(x,y2,':r','LineWidth',2);
semilogy(x,y3,'--','LineWidth',2);

xlim([0,5])
xlabel('x','FontSize',16,'FontWeight','Bold');
ylabel('PDF(x)','FontSize',16,'FontWeight','Bold');
title('Weibull densities','FontSize',16,'FontWeight','Bold');
ylim([10e-4 10e-1]);

set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on
% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss06_04.pdf
%print -painters -dpng -r600 STFloss06_04.png