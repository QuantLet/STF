function [x,y]=tstabpdf(alpha,sigma,lambda,beta,xmax,n,mult)
%linia 50 jest kluczowa
%STABPDF Stable probability density function (pdf).
%  [X,Y] = STABPDF(ALPHA,SIGMA,BETA,XMAX,N) returns the stable pdf ((X,Y) pairs) 
%  with characteristic exponent, ALPHA, scale, SIGMA, and skewness, BETA,
%  at 2^N evenly spaced values in [-XMAX,XMAX]. 
%  Default values for XMAX and N are 20 and 10 respectively.
%
%  For .95 < ALPHA < 1.05 the stable pdf is calculated using the formula 
%  for ALPHA = 1, otherwise numerical errors creep in. 
%
%  Due to the nature of FFT, values away from the center may be underestimated. 
%  For this reason STABPDF calculates the stable pdf on the interval 
%  [-XMAX,XMAX]*2^MULT and then truncates it to the original interval. 
%  The default value of MULT is 2, however, for better accuracy use MULT>4.
%  The full syntax is [X,Y] = STABPDF(ALPHA,SIGMA,BETA,XMAX,N,MULT).

%	Copyright (c) 1996-2001 by The Hugo Steinhaus Center
%	1996.07.01 Rafal Weron
%	2001.02.13 Rafal Weron

if nargin < 7,
   mult = 2;
end

if nargin < 6, 
    n = 10;
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
% if abs(alpha-1) < .05,
%   yy = exp(-sigma*(abs(xx)).*(1+i*beta.*sign(xx).*(2/pi).*log((abs(xx)))));
% else
%   yy = exp(-sigma^alpha.*(abs(xx)).^alpha.*(1-i*beta.*sign(xx).*tan((pi.*alpha)/2)));
% end;

%wzór z Matacza
%yy=exp(-sigma^alpha/cos(pi*alpha/2)*((xx.^2+lambda^2).^(alpha/2).*cos(alpha*atan(abs(xx)/lambda))-lambda^alpha));

%Matacz + Koponen, ¿eby by³ skoœny
yy=exp(-sigma^alpha/cos(pi*alpha/2)*(-lambda^alpha+(xx.^2+lambda^2).^(alpha/2).*cos(alpha*atan(abs(xx)/lambda)).*...
    (1+i*beta*sign(xx).*tan(alpha*atan(abs(xx)/lambda)))));


%wzór z Koponena !!! nie wiadomo czemu, ale gorszy
%yy=exp(lambda^alpha/cos(pi*alpha/2)-pi/(alpha*GAMMA(alpha)*sin(pi*alpha))*...
%    (xx.^2+lambda^2).^(alpha/2).*cos(alpha*atan(abs(xx)/lambda)));

%wzór z Wikipedii
%w innym pliku

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
