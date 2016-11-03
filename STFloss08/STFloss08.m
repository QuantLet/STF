
% clear variables and close windows
clear all
close all
clc


% load the data
data = load('dfl.dat','-ascii');
data = data(:,5);
xaxis = sort(data);

sampleMef = samplemef (data, xaxis);
sampleMef = sampleMef / 1e6;

xaxis = xaxis / 1e6;
figure
plot(xaxis,sampleMef,'ok');
ylabel('e_n(x) (DKK million)','FontSize',16,'FontWeight','Bold');
xlabel('x (DKK million)','FontSize',16,'FontWeight','Bold');
xlim([0,18]);
set(gca,'LineWidth',1.6,'FontSize',16,'FontWeight','Bold');
set(gca,'Ytick',[0:2:12],'YTickLabel',{0,2,4,6,8,10,12},'FontSize',16,'FontWeight','Bold')
box on

% to save the plot in pdf or png please uncomment next 2 lines:
%print -painters -dpdf -r600 STFloss08.pdf
%print -painters -dpng -r600 STFloss08.png
