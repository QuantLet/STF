
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **STF2vswapstrike** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of QuantLet : STF2vswapstrike

Published in : Statistical Tools for Finance and Insurance

Description : 'Calculates the strike of a variance of a given maturity using the potrfolio of
options of a given strike range.'

Keywords : variance, option, portfolio, maturity, strike

See also : STF2daxvswap, STF2dollargamma2D, STF2dollargamma3D, STF2logpayoff

Author : Elena Silyakova

Submitted : Mon, December 05 2011 by Dedy Dwi Prastyo

Usage : STF2vswapstrike

Input : 's1- the lowest strike of options" portfolio s2- the highest strike of options" portfolio
tau- maturity of the swap n- number of options used for replication'

Output : Strike value

```


### MATLAB Code:
```matlab
clear all
% user input parameters
disp('Please input the lowest and the highest strike of options portfolio as: [10,200]') ;
disp(' ') ;
para=input('Options strike range [lower bound, upper bound]=');
while length(para) < 2
  disp('Not enough input arguments. Please input in 1*2 vector form like [10,200] or [10,200]');
  para=input('Options strike range [lower bound, upper bound]=');
end
s1=para(1);
s2=para(2);

  disp(' ') ;
disp('Please imput the implied volatility, i.e. 0.25 for 25%') ;
disp(' ') ;
para=input('Implied volatility =');
v=para;

disp(' ') ;
disp('Please input the maturity of a swap (e.g. 1 year = 1, 3 month = 0.25)') ;
disp(' ') ;
para=input('Maturity of a swap =');
tau=para;

r=0;                        % interest rate
s=(s2+s1)/2;                % defines the underlying price S
kp=s1:10:(s-1);               % defines the range of puts
kc=(s+1):10:s2;               % defines the range of calls


callputspot = blsprice(s, s, r, tau, v); %BS value of the call(put) with strike equal to spot price


for  i=1:length(kc);   
ivc(i)=kc(i)*v/s;
d1=(log(s/kc(i))+(r+v^2/2)*tau)/(ivc(i)*tau^0.5);
d2=d1-ivc(i)*tau^0.5;

call(i) = s*normcdf(d1)-kc(i)*exp(-r*tau)*normcdf(d2);

i=i+1;    
end

for  i=1:length(kp)   
ivp(i)=kp(i)*v/s;

d1=(log(s/kp(i))+(r+ivp(i)^2/2)*tau)/(ivp(i)*tau^0.5);
d2=d1-ivp(i)*tau^0.5;
put(i) = kp(i)*exp(-r*tau)*normcdf(-d2)-s*normcdf(-d1);
i=i+1;
end    

fp=log(s./kp)+kp./s-1;
fc=log(s./kc)+kc./s-1;

wc(1)=(fc(1)-0)/(kc(1)-s); %weight of fist call(put) (with strike eqial to spot price)

wcum=wc(1);

for  j=2:length(kc);   
wc(j)=(fc(j)-fc(j-1))/(kc(j)-kc(j-1))-wcum;
wcum=wcum+wc(j);
j=j+1;
end

j=numel(kp);
wcum=(0-fp(j))/(kp(j)-s);

while kp(j)> s1   
wp(j)=(fp(j-1)-fp(j))/(kp(j)-kp(j-1))-wcum;
wcum=wcum+wp(j);
j=j-1;
end

wp = wp(2:length(wp));

put = put(2:length(put));
call = call(1:(length(call)-1));
opt=cat(2,put,callputspot,call);
w=cat(2,wp,wc);

strike=((2/tau)*sum(w.*opt))^0.5

clear x wcum  v th tau  r para opt n  j i fp f discr d2 d1 wp wc kp kc  fc




```
