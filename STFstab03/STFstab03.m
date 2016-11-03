% clear variables and close windows
clear all
close all
clc


cmd = [1 3];


if ismember(3,cmd),
    f = figure(3);
    subplot(2,2,1)
    [x1,y1] = stabpdf_fft(2,1,0,0,6);
    [x2,y2] = stabpdf_fft(1.7,1,0,0,20);
    x3 = linspace(-20,20,1000);
    y3 = tstabpdf3(x3,1.7,1,.2,0);
    semilogy(x1,y1,'b-','linewidth',2)
    hold on
    semilogy(x3,y3,'-.','color',[0 .5 0],'linewidth',1)
    semilogy(x2,y2,'r--','linewidth',1)
    hold off
    xlabel('x')
    ylabel('PDF(x)')
    set(gca,'xlim',[-10 10],'ylim',[1e-4 1]);

    subplot(2,2,2)
    loglog(x1,y1,'b-','linewidth',2)
    hold on
    loglog(x3,y3,'-.','color',[0 .5 0],'linewidth',1)
    loglog(x2,y2,'r--','linewidth',1)
    hold off
    set(gca,'xlim',[8e-1,3e1],'ylim',[1e-4 .5]);
    xlabel('x')
    ylabel('PDF(x)')
    legend('Gaussian','TSD(1.7,0.2)','1.7-Stable',1)

    subplot(2,2,3)
    y4 = tstabpdf3(x3,1.7,1,5,0);
    y5 = tstabpdf3(x3,1.7,1,.5,0);
    y6 = tstabpdf3(x3,1.7,1,.01,0);
    semilogy(x3,y4,'b-','linewidth',2)
    hold on
    semilogy(x3,y5,'k-','linewidth',1)
    semilogy(x3,y3,'-.','color',[0 .5 0],'linewidth',1)
    semilogy(x3,y6,'r--','linewidth',1)
    hold off
    xlabel('x')
    ylabel('PDF(x)')
    set(gca,'xlim',[-10 10],'ylim',[1e-4 1]);

    subplot(2,2,4)
    loglog(x3,y4,'b-','linewidth',2)
    hold on
    loglog(x3,y5,'k-','linewidth',1)
    loglog(x3,y3,'-.','color',[0 .5 0],'linewidth',1)
    loglog(x3,y6,'r--','linewidth',1)
    hold off
    set(gca,'xlim',[8e-1,3e1],'ylim',[1e-4 .5]);
    xlabel('x')
    ylabel('PDF(x)')
    legend('TSD(lambda=5)','TSD(lambda=0.5)','TSD(lambda=0.2)','TSD(lambda=0.01)',1)
    
    print(f,'-dpsc2','STF2stab03.ps')
end