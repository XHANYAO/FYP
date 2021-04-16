clc
clear all
GENMAX=10000000;
STATION=20000;
ACTIONNO=2000;
ALPHA=0.1;
GAMMA=0.9;
EPSILON=0.3;
nn=15;
GOAL=55;
LEVEL=nn;
INPUTNO=20000;
HIDDENNO=2000;
OUTPUTNO=2000;
NNALPHA=3;
count=0;
qvalue=0.999*rand(nn,2000);
x=[735534.5  652180.2  800773.8  715923.5  921134.9  892048.3  743773.8  790930.5  705540.3  869165.0  715923.5  921134.9  892048.3  743773.8  787157.0 984575.0 702192.0 768663.0 657394.0
   508229664 430814447 608749645 600559994 561337651 453478757 598871791 466361361 527207144 481861654 600559994 561337651 453478757 598871791 461197222 357248933 491403999 397368822 786297177
   3.052794  2.430809  2.423213  2.735298  2.517545  3.368810  2.838018  2.814051  2.229461  2.322614  2.735298  2.517545  3.368810  2.838018  1.802526 2.704607 3.348886 3.30791 2.217739
   2.925414  6.147478  3.530716  1.833023  2.653565  7.568968  5.010552  35.71505  4.165090  2.917472  1.833023  2.653565  7.568968  5.010552  3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565];

n(:,1)=(wgn(1,1000000,0));%信道噪声
hy=[n n n n n n n n n n n n n n n];%按理来说每个信道的噪声是不同的
y=(hy);
for n=1:nn
       m=1;
        for fi=0.2*10^9:(0.8*1e+9)/999:1*10^9
            t(n,m)=x(2,n)/fi;
            if t(n,m)>x(3,n)
                t(n,m)=10000;
            end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for fi=0.2*10^9:(0.8*1e+9)/999:1*10^9
            e(n,m)=x(2,n)*fi^2*10^(-26);
            m=m+1;
        end
    
end
GL=0.2*t+0.8*e;
for n=1:nn
gl(n)=min(GL(n,:));
end
cl=sum(gl)
cost=cl;
for n=1:nn
[rowl(n) coll(n)]=find(GL(n,:)==min(GL(n,:)));

TL(n)=t(n,coll(n));

EL(n)=e(n,coll(n));
end
%服务器时间和能耗（无噪声）
for n=1:nn
       m=1;
        for pi=0.001:0.1985/999:10^-0.7
           
            ri=0.2*10^6*log2(1+pi*x(4,n)^2);%上行传输速率
            
            tc(n,m)=x(1,n)/ri+x(2,n)/(4*10^9);
           if tc(n,m)>x(3,n)
                tc(n,m)=10000;
           end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for pi=0.001:0.1985/999:10^-0.7
            ri=0.2*10^6*log2(1+pi*x(4,n)^2);
            ec(n,m)=pi*x(1,n)/ri;
            m=m+1;
        end
    
end
GC=0.2*tc+0.8*ec;
for n=1:nn
gc(n)=min(GC(n,:));
end
for n=1:nn
[rowc(n) colc(n)]=find(GC(n,:)==min(GC(n,:)));

TC(n)=tc(n,colc(n));

EC(n)=ec(n,colc(n));
end
g=[gl;gc];
for n=1:nn
G(n)=min(g(:,n));
end
C=sum(G) %无噪
%initialize Q value，hide h，output
for n=1:nn
    if g(1,n)<g(2,n)
    qvalue(n,coll(n))=1;
    else
    qvalue(n,colc(n)+1000)=1;
    end
end
for i=1:HIDDENNO
    for j=1:INPUTNO+1
        wh(i,j)=1-2*rand();
    end
end
i=1;j=1;
for i=1:OUTPUTNO
    for j=1:HIDDENNO+1
        wo(i,j)=1-2*rand();
    end
end
i=1;j=1;
for i=1:GENMAX
    s=0;%initial state
    m=1;
    
    for t=1:LEVEL
        ran=rand();
        if ran<EPSILON%choose action 
        a= randi(2000);
            if a<=1000
               fi=(0.8*rand()+0.2)*(1e+9);
               ti=x(2,t)/fi;
               ei=x(2,t)*fi^2*10^(-26);
               c=0.2*ti+0.8*ei;
            else
               pi=0.1995*rand();
               for l=1:100
               r(l)=0.2*10^6*log2(1+pi*(x(4,t)+y(l,t))^2);
               end
               ri=mean(r);
               ti=x(1,t)/ri+x(2,t)/(4*10^9);
               ei=pi*x(1,t)/ri;
               c=0.2*ti+0.8*ei;
            end
        else
            [rowl coll]=find(qvalue(t,:)==max(qvalue(t,:)));
            a=coll;
            if a<=1000
               fi=0.2*1e+9+coll*(0.8*1e+9)/999;
               ti=x(2,t)/fi;
               ei=x(2,t)*fi^2*10^(-26);
               c=0.2*ti+0.8*ei;
            else
                pi=0.001+(coll-1000)*0.1985/999;
                for l=1:100
                r(l)=0.2*10^6*log2(1+pi*(x(4,t)+y(l,t))^2);
                end
                ri=mean(r);
               ti=x(1,t)/ri+x(2,t)/(4*10^9);
               ei=pi*x(1,t)/ri;
               c=0.2*ti+0.8*ei;
            end
        end
        snext=s+c;
       %[rownext colnext]=find(qvalue(t+1,:)==max(qvalue(t+1,:)));
        %anext=colnext;
        REWARD=10*(1-c/cl);%cl:local total cost
        erro=(c-G(t))^2;
        if erro<=0.0001%update Q
        qv=qvalue(t,a)+ALPHA*(REWARD-qvalue(t,a));
        else 
        qv=qvalue(t,a)+ALPHA*(REWARD/3-qvalue(t,a));
        end
    
    qvalue(t,a)=qv;
    s=snext;
    end
    if s<cost
        cost=s;
    end
    ccc=(C-cost)^2;
    if ccc<0.002
        break
    end
    ans=qvalue;%output Q
end
            
ttcc=[5.8034 6.6017 7.3395 7.7909 8.3374 8.9423] ;
