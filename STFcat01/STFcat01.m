clear all
close all
clc
d= load('ncl.dat')
plot(1990+d(:,2),d(:,3)/1e9,'LineWidth',2)
xlabel('Years','fontsize',12)
ylabel('Adjusted PCS catastrophe claims (USD billion)','fontsize',12)
