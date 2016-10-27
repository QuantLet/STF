% clear variables and close windows
clear all
close all
clc


cmd = [1 4];

if ismember(4,cmd),
    f = figure(4)
    x = linspace(-10,10,1000);
    % Make all distributions have the same variance
    [evh,varh] = ghypstat(1,3,0,.444,0);
    [evn,varn] = ghypstat(-0.5,1.5,0,.5,0);
    alpha = 1.7; lambda = 0.5; sigma = .391;
    vartsd = (alpha*(1-alpha)/cos(pi*alpha/2))*sigma^alpha*lambda^(alpha-2);
    [varh varn vartsd]  
    y1 = hyppdf(x',3,0,.444,0);
    y2 = nigpdf(x',1.5,0,.5,0);
    y3 = normpdf(x,0,varn^.5);
    y4 = tstabpdf3(x,alpha,sigma,lambda,0);
        
    subplot(1,2,1)
    plot(x,y1,'b-','linewidth',2)
    hold on
    plot(x,y4,'k-','linewidth',1)
    plot(x,y2,'-.','color',[0 .5 0],'linewidth',1)
    plot(x,y3,'r--','linewidth',1)
    hold off
    xlabel('x')
    ylabel('PDF(x)')
    set(gca,'xlim',[0,2]);
    legend('Hyperbolic','TSD(1.7,0.5)','NIG','Gaussian')

    subplot(1,2,2)
    semilogy(x,y1,'b-','linewidth',2)
    hold on
    semilogy(x,y4,'k-','linewidth',1)
    semilogy(x,y2,'-.','color',[0 .5 0],'linewidth',1)
    semilogy(x,y3,'r--','linewidth',1)
    hold off
    set(gca,'xlim',[-10,10],'ylim',[10e-9 1])
    xlabel('x')
    ylabel('PDF(x)')

    print(f,'-dpsc2','STF2stab04.ps')
end