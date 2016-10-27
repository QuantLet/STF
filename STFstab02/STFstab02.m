% clear variables and close windows
clear all
close all
clc


cmd = [1 2];


if ismember(2,cmd),
    f = figure(2);
    subplot(1,2,1)
    [x1,y1] = stabpdf_fft(2,1,0,0,6);
    [x2,y2] = stabpdf_fft(1.9,1,0,0,10);
    [x3,y3] = stabpdf_fft(1.5,1,0,0,10);
    [x4,y4] = stabpdf_fft(0.5,1,0,0,10);
    semilogy(x1,y1,'b-','linewidth',2)
    hold on
    semilogy(x2,y2,'-.','color',[0 .5 0],'linewidth',1)
    semilogy(x3,y3,'r--','linewidth',1)
    semilogy(x4,y4,'k-','linewidth',1)
    hold off
    xlabel('x')
    ylabel('PDF(x)')
    legend('alpha=2','alpha=1.9','alpha=1.5','alpha=0.5','Location','South')
    set(gca,'ylim',[1e-4 .8]);

    subplot(1,2,2)
    [x1,y1] = stabpdf_fft(1.2,2,0,0,5);
    [x2,y2] = stabpdf_fft(1.2,2,-1,0,5);
    [x3,y3] = stabpdf_fft(1.2,2,0.5,0,5);
    [x4,y4] = stabpdf_fft(1.2,2,1,0,5);
    plot(x1,y1,'b-','linewidth',2)
    hold on
    plot(x2,y2,'-.','color',[0 .5 0],'linewidth',1)
    plot(x3,y3,'r--','linewidth',1)
    plot(x4,y4,'k-','linewidth',1)
    hold off
    xlabel('x')
    ylabel('PDF(x)')
    legend('beta=0','beta=-1','beta=0.5','beta=1','Location','South')
    set(gca,'xtick',-5:2.5:5,'ylim',[0 .16])
    
    print(f,'-dpsc2','STF2stab02.ps')
end