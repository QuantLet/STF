% clear variables and close windows
clear all
close all
clc


xaxis = (1:2000) / 100;

% log
p1 = mef (0.001, 1, -1, xaxis,1);
% gamma 1
p2 = mef (0.5, 1/1.63510, -1, xaxis,2);
% gamma 2
p3 = mef (1.5, 1/1.63510, -1, xaxis,2);
% mix
p4 = mef (1/4, 1/0.8,0.1, xaxis,6);

plot(xaxis,p1,'-k','LineWidth',2);
hold on
plot(xaxis,p2,'--','LineWidth',2);
plot(xaxis,p3,'-.g','LineWidth',2);
plot(xaxis,p4,':r','LineWidth',2);
ylabel('e(x)','FontSize',16,'FontWeight','Bold');
xlabel('x','FontSize',16,'FontWeight','Bold');
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss07_01.pdf
%print -painters -dpng -r600 STFloss07_01.png

% MEF for Pareto, Burr, 2 x Weibull

% Weibull
p1 = mef (1.1, 0.5, -1, xaxis,3);
% Weibull
p2 = mef (0.9, 0.5, -1, xaxis,3);
% Pareto
p3 = mef (5 , 3.03, -1, xaxis,4);
%Burr
p4 = mef (2.39, 3.03,3, xaxis,5);

figure
plot(xaxis,p1,'--g','LineWidth',2);
hold on
plot(xaxis,p2,'-.b','LineWidth',2);
plot(xaxis,p3,'-k','LineWidth',2);
plot(xaxis,p4,':r','LineWidth',2);
ylabel('e(x)','FontSize',16,'FontWeight','Bold');
xlabel('x','FontSize',16,'FontWeight','Bold');
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss07_02.pdf
%print -painters -dpng -r600 STFloss07_02.png
