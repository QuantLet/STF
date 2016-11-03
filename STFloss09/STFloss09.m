
% clear variables and close windows
clear all
close all
clc


m = load('dfl.dat');
v = m(find(m(:,5)~=0),5:6);

n = getQnumber(v(:,2));

%  Qnumbers
t = (1:length(n))/4;

plot(t,n, 'LineWidth',2);
xlabel('Time (years)','FontSize',16,'FontWeight','Bold');
ylabel('Number of events','FontSize',16,'FontWeight','Bold'); 
xlim([0,23]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss09_01.pdf
%print -painters -dpng -r600 STFloss09_01.png

%  Mean value function

ncum = cumsum(n);
figure
plot(t,ncum,'--', 'LineWidth',2);
hold on
plot(t,(100.2394*t),'k', 'LineWidth',2);
plot(t,17.9937*t+3.5759*t.^2,'-.r', 'LineWidth',2)
% ylim([0,700])
xlim([0,23]);
xlabel('Time (years)','FontSize',16,'FontWeight','Bold');
ylabel('Aggregate number of losses / Mean value function','FontSize',16,'FontWeight','Bold');

set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss09_02.pdf
%print -painters -dpng -r600 STFloss09_02.png

