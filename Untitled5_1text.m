clc
clear all
%Ìõ¼þ
x=[735534.5  652180.2  800773.8  715923.5  921134.9  892048.3  743773.8  790930.5  705540.3  869165.0  715923.5  921134.9  892048.3  743773.8  787157.0 984575.0 702192.0 768663.0 657394.0
   508229664 430814447 608749645 600559994 561337651 453478757 598871791 466361361 527207144 481861654 600559994 561337651 453478757 598871791 461197222 357248933 491403999 397368822 786297177
   3.052794  2.430809  2.423213  2.735298  2.517545  3.368810  2.838018  2.814051  2.229461  2.322614  2.735298  2.517545  3.368810  2.838018  1.802526 2.704607 3.348886 3.30791 2.217739
   2.925414  6.147478  3.530716  1.833023  2.653565  7.568968  5.010552  35.71505  4.165090  2.917472  1.833023  2.653565  7.568968  5.010552  3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565];

n(:,1)=(wgn(1,1000000,0));
hy=[n n n n n n n n n n n n n n n];
y=hy;
t=5;
k=10^(-26);
for nn=1:15
     no(:,nn)=(wgn(1,1000000,0));%channel error
    y=no;
%local time and energy
for n=1:nn
       m=1;
        for fi=0.2*10^9:0.04*10^8:1*10^9
            tl(n,m)=x(2,n)/fi;%local compute time
            if tl(n,m)>x(3,n)
                tl(n,m)=100;
            end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for fi=0.2*10^9:0.04*10^8:1*10^9
            el(n,m)=x(2,n)*fi^2*10^(-26);%local energy calculate
            m=m+1;
        end
    
end
GL=0.2*tl+0.8*el;
for n=1:nn
gl(n)=min(GL(n,:));
end
%local minimum total cost
for n=1:nn
[rowl(n) coll(n)]=find(GL(n,:)==min(GL(n,:)));
end
for n=1:nn
TL(n)=tl(n,coll(n));
end
for n=1:nn
EL(n)=el(n,coll(n));
end
%MEC time and energy
for n=1:nn
       m=1;
        for pi=0.001:0.00099:0.2
             for l=1:500
            rr(l)=0.2*10^6*log2(1+pi*(x(4,n)+y(l,n))^2);
            end
            ri=mean(rr);%upload transform speed
            tc(n,m)=x(1,n)/ri+x(2,n)/(4*10^9);
           if tc(n,m)>x(3,n)
                tc(n,m)=100;
           end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for pi=0.001:0.00099:0.2
             for l=1:500
            rr(l)=0.2*10^6*log2(1+pi*(x(4,n)+y(l,n))^2);
            end
            ri=mean(rr);
            ec(n,m)=pi*x(1,n)/ri;
            m=m+1;
        end
    
end

GC=0.2*tc+0.8*ec;
for n=1:nn
gc(n)=min(GC(n,:));
end
%{
for n=1:10
gcc(n)=mean(gc(:,n));
end
%}
%MEC minimum total cost
for n=1:nn
[rowc(n) colc(n)]=find(GC(n,:)==min(GC(n,:)));
end
for n=1:nn
TC(n)=tc(n,colc(n));
end
for n=1:nn
EC(n)=ec(n,colc(n));
end

%total cost 
g=[gl;gc];
for n=1:nn
G(n)=min(g(:,n));
end
C(nn)=sum(G)

%ESO
%{
for n=1:10
    for m=1:201
EI(n,m)=min([el(n,m);ec(n,m)]);
    end
end
for n=1:10
 for m=1:201
TI(n,m)=min([tl(n,m);tc(n,m)]);
    end
end
%}
for n=1:nn
r(n,1)=1;
u(n,1)=1;
b=0.1;
ao=1e+10;
gr(n,1)=1;
gu(n,1)=1;

eexc(n,1)=0;
etxc(n,1)=0;
st(n,1)=0;
sp(n,1)=0;
for m=1:100000

     
fx=(u(n,1)/2/k/r(n,1))^(1/3);
if fx<0.2*10^9
    fxI(n,1)=0.2*10^9;
else if fx>1*10^9
         fxI(n,1)=1*10^9;
    else fxI(n,1)=real(fx);
    end
end
p=1;
op=0;
o=inf;
for pi=0.001:0.0001:0.2
    ri=0.2*10^6*log2(1+pi*(x(4,n)+hy(m,n))^2);
    v(p)= x(1,n)*(r(n,1)*pi+u(n,1))/ri; 
   
       if v(p)<o
        o=v(p);
        op=pi;
       end
   
    p=p+1;
end
%{
o=gradient(v);
[rowp,colp]=find(o==max(o));
op=0,001+0.0001*(colp-1);
%}
pmin=(2^(x(1,n)/(0.2*10^6)/(x(3,n)-x(2,n)/(4*10^9)))-1)/((x(4,n)+hy(m,n))^2);
if pmin>0.2;
    pmin=0.2;
end
if op<pmin
    pxi(n,1)=pmin;
else if op>0.2
         pxi(n,1)=0.2;
    else pxi(n,1)=op;
    end
end

txl(n,1)=x(2,n)/fxI(n,1);
exl(n,1)=x(2,n)*fxI(n,1)^2*10^(-26);
GxL(n,1)=0.2*txl(n,1)+0.8*exl(n,1);
for l=1:500
     rx(l)=0.2*10^6*log2(1+pxi(n,1)*(x(4,n)+hy(l,n))^2);
end    
rxi(n,1)=mean(rx);
txc(n,1)=x(1,n)/rxi(n,1)+x(2,n)/(4*10^9);
exc(n,1)=pxi(n,1)*x(1,n)/rxi(n,1);
GxC(n,1)=0.2*txc(n,1)+0.8*exc(n,1);
pmax=0.2;
for l=1:500
rma(l)=0.2*10^6*log2(1+pmax*(x(4,n)+hy(l,n))^2);
end
rmax=mean(rma);
tcmin(n)=x(1,n)/rmax+x(2,n)/(4*10^9);
if GxL(n,1)<GxC(n,1)
    s(n,1)=0;
    GG(n,1)=GxL(n,1);
    for TIL=txl(n,1):0.001:x(3,n)
        EIL=k*x(2,n)^3/TIL;
     a=-r(n,1)*EIL-u(n,1)*TIL+0.8*EIL+0.2*TIL;
     if a<ao 
          ao=a;
          EIO(n,1)=EIL;
          TIO(n,1)=TIL;
     end
    end
else
    s(n,1)=1;
    GG(n,1)=GxC(n,1);
for TI=tcmin(n):0.001:txc(n,1)
        EI=(TI-x(2,n)/(4*10^9))*(2^(x(1,n)/(0.2*10^6)/(TI-x(2,n)/(4*10^9)))-1)/((x(4,n)+hy(m,n)^2));
       a=-r(n,1)*EI-u(n,1)*TI+0.8*EI+0.2*TI;
     if a<ao
          ao=a;    
          EIO(n,1)=EI;
          TIO(n,1)=TI;
     end
end
end
if s(n,1)==1
eexc(n,1)=(eexc(n,1)+exc(n,1));
etxc(n,1)=(etxc(n,1)+txc(n,1));
end
gr(n,1)=EIO(n,1)-s(n,1)*exc(n,1)-(1-s(n,1))*exl(n,1);
gu(n,1)=TIO(n,1)-s(n,1)*txc(n,1)-(1-s(n,1))*txl(n,1);
%{
sp(n,1)=gr(n,1)*r(n,1)+gu(n,1)*u(n,1);
st(n,1)=(st(n,1)+sp(n,1));
stp=st(n,1)/m;
if (-0.01<=st(n,1)/m)&&(st(n,1)/m<=0.01)
    break
end
%}
r(n,1)=abs(r(n,1)-b*gr(n,1));
u(n,1)=abs(u(n,1)-b*gu(n,1));
GXI(n,1)=0.8*(EIO(n,1))+0.2*(TIO(n,1));
finalc=G(n);
realc=GG(n,1);
stop=abs(finalc-realc);
if stop<0.001
    break
end
end
ANS=sum(GG);
end
TOLC(nn)=ANS;
end
plot(TOLC)