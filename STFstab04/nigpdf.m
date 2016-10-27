function y=nigpdf(x,a,b,d,m)
if a^2-b^2>0 & a>0 & a<100
    g=sqrt(a^2-b^2);
    y=a*d*besselk(1,a*sqrt(d^2+(x-m).^2))./(pi*sqrt(d^2+(x-m).^2)).*exp(d*g+b*(x-m));
else
    y=NaN(size(x));
end;