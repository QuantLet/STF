function retval = umlp(x,y)
% unidirextional minimum length path algorithm
% x - the distance matrix
% y - root of the chain (the number of column)
%
% The rsult is presentes as a set of links between nodes
%
% Author: Janusz Miśkiewicz, email: jamis@ift.uni.wroc.pl
[n,m]=size(x);
net=zeros(n-1,3);
onnet=zeros(n,1);
licz=0;
onnet(y)=1;
maxx=10*max(max(x));
smax=maxx*eye(size(x));
x=x+smax;
while (licz<n-1)
  minx=min(x(y,:));
  it=find(x(y,:)==minx);
    if (length(it) > 1)
      tmp=1;
      while ((onnet(it(tmp))==1) && (tmp < length(it)))
	   tmp=tmp+1;
      end;
      if (onnet(it(tmp))==0)
        ii=it(tmp);
        it=[];
        it=ii;
        tmp=0;
      else
        ii=it(1);
        it=[];
        it=ii;
      end;
    end;
    if (onnet(it)==0 )
      licz=licz+1;
      net(licz,1)=y;
      net(licz,2)=it;
      net(licz,3)=x(it,y);
      onnet(it)=1;
      x(it,y)=maxx;
      x(y,it)=maxx;
      y=it;
    end;
    
    if ((onnet(it)==1) && (onnet(y)==1))
      x(it,y)=maxx;
      x(y,it)=maxx;
    end;
end;
retval=net;
end