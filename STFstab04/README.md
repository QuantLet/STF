
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STFstab04** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of QuantLet : STFstab04

Published in : Statistical Tools for Finance and Insurance

Description : 'Plots densities and log-densities of symmetric hyperbolic, TSD, NIG and Gaussian
distributions having the same variance. Requires ghypstat.m (requires besselk.m from Symbolic Math
toolbox), hyppdf.m, nigpdf.m, tstabpdf3.m, and tstabpdf.m functions.'

Keywords : 'pdf, density, hyperbolic, stable distribution, gaussian, graphical representation,
visualization'

See also : ghypstat, hyppdf, nigpdf, tstabpdf, tstabpdf3

Author : Rafal Weron

Submitted : Tue, September 18 2012 by Dedy Dwi Prastyo

Example : ghypstat.m, hyppdf.m, nigpdf.m, tstabpdf3.m, tstabpdf.m produce these plots.

```

![Picture1](plot.png)


### MATLAB Code:
```matlab
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
```
