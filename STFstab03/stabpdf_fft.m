function [x,y]=stabpdf_fft(alpha,sigma,beta,mu,xmax,n,S0,mult)
%STABPDF_FFT Stable probability density function (pdf) via FFT.
%   [X,Y] = STABPDF_FFT(ALPHA,SIGMA,BETA,MU,XMAX,N) returns the stable pdf 
%   ((X,Y) pairs) with characteristic exponent, ALPHA, scale, SIGMA, 
%   skewness, BETA, and location MU at 2^N evenly spaced values in 
%   [-XMAX,XMAX] in the S0 parametrization. Default values for XMAX and N 
%   are 20 and 12 respectively.
%   [X,Y] = STABPDF_FFT(ALPHA,SIGMA,BETA,MU,XMAX,N,0) returns the stable 
%   pdf in the S parametrization.
%
%   For .999 < ALPHA < 1.001 the stable pdf is calculated using the formula 
%   for ALPHA = 1, otherwise numerical errors creep in. 
%
%   Due to the nature of FFT, values away from the center may be underestimated. 
%   For this reason STABPDF_FFT calculates the stable pdf on the interval 
%   [-XMAX,XMAX]*2^MULT and then truncates it to the original interval. 
%   The default value of MULT is 4, however, for better accuracy use MULT>4.
%   The full syntax is [X,Y] = STABPDF_FFT(ALPHA,SIGMA,BETA,XMAX,N,PARAM,MULT).
%
%   Reference(s):
%	[1] R.Weron (2004) "Computationally intensive Value at Risk 
%   calculations", in "Handbook of Computational Statistics: Concepts and 
%   Methods", eds. J.E. Gentle, W. Haerdle, Y. Mori, Springer, Berlin, 
%   911-950. 

%	Copyright (c) 1996-2009 by The Hugo Steinhaus Center
%	1996.07.01 Rafal Weron
%	2001.02.13 Rafal Weron
%   2009.03.27 Rafal Weron

if nargin < 8,
   mult = 4;
end

if nargin < 7,
   S0 = 1;
end

if nargin < 6, 
    n = 12;
end

if nargin < 5;
    xmax = 20;
end

% calculate pdf on a twice larger interval
n = n+mult;
xmax = xmax*(2^mult);

M = 2^n;
R = pi/xmax;
dt = 1/(R*M);
  
xx = (-2^(n-1)+.5:(2^(n-1)-.5))/(2^n*dt);

% stable characteristic function
if abs(alpha-1) < .001,
    if S0, 
        mu = mu - beta*sigma*2/pi*log(sigma);
    end
    yy = exp(-sigma*(abs(xx)).*(1+i*beta.*sign(xx).*(2/pi).*log((abs(xx)))) + i*mu*xx);
else
    if S0,
        mu = mu - beta*sigma*tan(0.5*pi*alpha);
    end
    yy = exp(-sigma^alpha.*(abs(xx)).^alpha.*(1-i*beta.*sign(xx).*tan((pi.*alpha)/2)) + i*mu*xx);
end;
  
% FFT
yy1 = [yy((2^(n-1)+1):2^n), yy(1:2^(n-1))];
z = real( fft(yy1) )/(2*pi)*R;

% stable density
x = (2*pi)*((0:1:(M-1))/(M*R)-1/(2*R));
y = [z((2^(n-1)+1):2^n), z(1:2^(n-1))];   

% shrink to the original interval
T = find((x<=xmax/(2^mult)) & (x>=-xmax/(2^mult)));
x = x(T);
y = y(T);
