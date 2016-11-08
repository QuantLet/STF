function retval = mst(x)
% Algorithm generates minimum spanning tree
% The rsult is presentes as a set of links between nodes
%
% Author: Janusz Miśkiewicz, email: jamis@ift.uni.wroc.pl
[n,m]=size(x);
x=triu(x,1);
net=zeros(n-1,3);
onnet=zeros(n,1);
klaster=zeros(n,1);
klast=0;
licz=0;
%check if the matrics is symmetric and positive
maxx=max(max(x));
smax=10*abs(maxx);
x(x==0)=smax;
while (licz<n-1)
  minx=min(min(x));
  [i,j]=find(x<=minx);
  if (length(i) > 0) 
    ii=i(1);
    jj=j(1);
    i=[];
    j=[];
    i=ii;
    j=jj;
  end;
  if (onnet(i) ==0 && onnet(j) ==0)
    licz=licz+1;
    net(licz,1)=i;
    net(licz,2)=j;
    klast=klast+1;
    klaster(i)=klast;
    klaster(j)=klast;
    net(licz,3)=min(x(i,j),x(j,i));
    onnet(i)=1;
    onnet(j)=1;
    x(i,j)=smax;
    x(j,i)=smax;
  elseif (onnet(i)==0 && onnet(j)==1)  
    licz=licz+1;
    net(licz,1)=i;
    net(licz,2)=j;
    net(licz,3)=min(x(i,j),x(j,i));
    onnet(i)=1;
    klaster(i)=klaster(j);
    x(i,j)=smax;
    x(j,i)=smax;
  elseif (onnet(i) ==1 && onnet(j) ==0)  
    licz=licz+1;
    net(licz,1)=i;
    net(licz,2)=j;
    net(licz,3)=min(x(i,j),x(j,i));
    onnet(j)=1;
    klaster(j)=klaster(i);
    x(i,j)=smax;
    x(j,i)=smax;
  elseif (onnet(i) ==1 && onnet(j) ==1 && klaster(i)==klaster(j))  
    x(i,j)=smax;
    x(j,i)=smax;
  elseif  (onnet(i) ==1 && onnet(j) ==1 && klaster(i)~=klaster(j))
    licz=licz+1;
    net(licz,1)=i;
    net(licz,2)=j;
    net(licz,3)=min(x(i,j),x(j,i));
    klaster(klaster==klaster(i))=klaster(j);
  end;
end;
retval=net;
end