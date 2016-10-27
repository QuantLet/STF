% clear variables and close windows
clear all
close all
clc

cmd = [1 1];

dane = load('Dane_copula.txt');



if ismember(1,cmd),
    %d - DJIA
    %x - zwroty z d
    %xx - linspace na zwrotach d
    d = dane(:,7); % DJIA data
    d = d(~isnan(d));
    x = log(d(2:end)./d(1:end-1));
    xx = linspace(min(x),max(x),1000);
    f = figure(1);
    subplot(1,2,1)
    plot(x)
    xlabel('Days (2000.01.03-2009.12.31)')
    ylabel('DJIA Returns')
    set(gca,'xlim',[0,length(x)+1],'ylim',[-.12 .12])

    subplot(1,2,2)
    [xemp,yemp] = empcdf(x);
    [mu,sig] = normfit(x);
    loglog(-xemp,yemp,'k.')
    hold on
    loglog(-xx,normcdf(xx,mu,sig),'r--','linewidth',1)
    hold off
    xlabel('-x')
    ylabel('CDF(x)')
    set(gca,'xlim',[1e-3 1e-1],'ylim',[1e-4 1])
    legend('DJIA returns','Gaussian fit',3)
    
    print(f,'-dpsc2','STF2stab01.ps')
end