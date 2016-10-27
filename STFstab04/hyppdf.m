function y=hyppdf(x,a,b,d,m)
if a^2-b^2>0 & a>0
    g=sqrt(a^2-b^2);
    y=g/(2*a*d*besselk(1,d*g))*exp(-a*sqrt(d^2+(x-m).^2)+b*(x-m));
else
    y=NaN(size(x));
end;