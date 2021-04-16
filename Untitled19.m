clc
clear all
GENMAX=200000;
STATION=21;
ACTIONNO=2;
ALPHA=0.1;
GAMMA=0.9;
EPSILON=0.7;
nn=15;
GOAL=55;
LEVEL=nn;
INPUTNO=2;
HIDDENNO=2;
OUTPUTNO=1;
NNALPHA=3;
count=0;
qvalue=rand(nn,2);
x=[735534.5  652180.2  800773.8  715923.5  921134.9  892048.3  743773.8  790930.5  705540.3  869165.0  715923.5  921134.9  892048.3  743773.8  787157.0 984575.0 702192.0 768663.0 657394.0
   508229664 430814447 608749645 600559994 561337651 453478757 598871791 466361361 527207144 481861654 600559994 561337651 453478757 598871791 461197222 357248933 491403999 397368822 786297177
   3.052794  2.430809  2.423213  2.735298  2.517545  3.368810  2.838018  2.814051  2.229461  2.322614  2.735298  2.517545  3.368810  2.838018  1.802526 2.704607 3.348886 3.30791 2.217739
   2.925414  6.147478  3.530716  1.833023  2.653565  7.568968  5.010552  35.71505  4.165090  2.917472  1.833023  2.653565  7.568968  5.010552  3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565];

n(:,1)=(wgn(1,1000000,0));%信道噪声
hy=[n n n n n n n n n n n n n n n];%按理来说每个信道的噪声是不同的
y=(hy);
for n=1:nn
       m=1;
        for fi=0.2*10^9:0.001*10^8:1*10^9
            t(n,m)=x(2,n)/fi;
            if t(n,m)>x(3,n)
                t(n,m)=inf;
            end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for fi=0.2*10^9:0.001*10^8:1*10^9
            el(n,m)=x(2,n)*fi^2*10^(-26);
            m=m+1;
        end
    
end
GL=0.2*t+0.8*el;
for n=1:nn
gl(n)=min(GL(n,:));
end
cl=sum(gl)
cost=cl;
%服务器时间和能耗（无噪声）
for n=1:nn
       m=1;
        for pi=0.001:0.0005:10^-0.7
            for l=1:500
            r(l)=0.2*10^6*log2(1+pi*(x(4,n)+y(l,n))^2);
            end
            ri=mean(r);
           
            tc(n,m)=x(1,n)/ri+x(2,n)/(4*10^9);
           if tc(n,m)>x(3,n)
                tc(n,m)=inf;
           end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for pi=0.001:0.0005:10^-0.7
            for l=1:500
            r(l)=0.2*10^6*log2(1+pi*(x(4,n)+y(l,n))^2);
            end
            ri=mean(r);
            ec(n,m)=pi*x(1,n)/ri;
            m=m+1;
        end
    
end
GC=0.2*tc+0.8*ec;
for n=1:nn
gc(n)=min(GC(n,:));
end
g=[gl;gc];
for n=1:nn%teacher data
if gl(n)<gc(n)
    e(n,1)=0;
    e(n,2)=1;
    e(n,3)=0;
else
    e(n,1)=1;
    e(n,2)=0;
    e(n,3)=1;
end
end
    
for n=1:nn
G(n)=min(g(:,n));
end
C=sum(G) %noise
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
for loop=1:GENMAX%Q value
    s=0;
    m=1;
    
    for t=1:LEVEL
        ran=rand();
        if ran<EPSILON%choose action 
        a= randi(2);
        crl=inf;
               for fi=0.2*10^9:0.001*10^8:1*10^9
                   til=x(2,t)/fi;
                   if til>x(3,t)
                      til=100000;
                    end
                   eil=x(2,t)*fi^2*10^(-26);
                   gil=0.2*til+0.8*eil;
                   
                   if gil<crl
                       crl=gil;
                   end
               end
                crc=inf;
               for pi=0.001:0.0005:10^-0.7
                 for l=1:500
                    rir(l)=0.2*10^6*log2(1+pi*(x(4,t)+y(l,t))^2);
                  end
                    ric=mean(rir);
                 
               tic=x(1,t)/ric+x(2,t)/(4*10^9);
               if tic>x(3,t)
                      tic=100000;
               end
               eic=pi*x(1,t)/ric;
               gic=0.2*tic+0.8*eic;
          
                   if gic<crc
                       crc=gc(t);
                   end
               end
           if a==1
               c=crl;
           else
               c=crc;
           end
           if crl<crc
               rc=crl;
           else
               rc=crc;
           end         
        else
            for i=1:HIDDENNO% 
                u(i)=0;
                for j=1:INPUTNO
                    u(i)=u(i)+e(t,j)*wh(i,j);
                end
                u(i)=u(i)-wh(i,INPUTNO+1);
                hi(i)=1/(1+exp(-u(i)));
            end
            for num=1:OUTPUTNO
            o(num,t)=0;
            for i=1:HIDDENNO
                o(num,t)=o(num,t)+hi(i)*wo(num,i);
            end
            o(num,t)=o(num,t)-wo(num,i);
            o(num,t)=1/(1+exp(-o(num,t)));
            end
            
            a=round(o(num,t)+1);%action choise
             cfl=inf;
             for fi=0.2*10^9:0.001*10^8:1*10^9
                   til=x(2,t)/fi;
                   if til>x(3,t)
                      til=100000;
                    end
                   eil=x(2,t)*fi^2*10^(-26);
                   gil=0.2*til+0.8*eil;
                  
                   if gil<cfl
                       cfl=gil;
                   end
             end
               cfc=inf;
                for pi=0.001:0.0005:10^-0.7
                 for l=1:500
                    rir(l)=0.2*10^6*log2(1+pi*(x(4,t)+y(l,t))^2);
                  end
                    ric=mean(rir);
               tic=x(1,t)/ric+x(2,t)/(4*10^9);
               if tic>x(3,t)
                      tic=100000;
               end
               eic=pi*x(1,t)/ric;
               gic=0.2*tic+0.8*eic;
                  
                   if gic<cfc
                       cfc=gc(t);
                   end
               end
            if a==1
               c=cfl;
           else
               c=cfc;
           end
           if cfl<cfc
               rc=cfl;
           else
               rc=cfc;
           end        
                
        end
        
        snext=s+c;
       %[rownext colnext]=find(qvalue(t+1,:)==max(qvalue(t+1,:)));
        %anext=colnext;
        REWARD=(1-c/cl);
        if c==rc
        qv=qvalue(t,a)+ALPHA*(3*REWARD-qvalue(t,a));
        else
        qv=qvalue(t,a)+ALPHA*(REWARD-qvalue(t,a));
        end
        
    qvalue(t,a)=qv;%upload Q 
       for i=1:HIDDENNO%neural network
           u(i)=0;
           for n=1:INPUTNO
           u(i)=u(i)+e(t,n)*wh(i,n);
           end
           u(i)=u(i)-wh(i,3);
           hi(i)=1/(1+exp(-u(i)));
       end
       for num=1:OUTPUTNO
       o(num,t)=0;
       for i=1:HIDDENNO
         o(num,t)=o(num,t)+hi(i)*wo(num,i);
       end
        o(num,t)=o(num,t)-wo(num,3);
        o(num,t)=1/(1+exp(-o(num,t)));
       end
       d=(e(t,INPUTNO+num)-o(num,t))*o(num,t)*(1-o(num,t));
       for m=1:HIDDENNO
           wo(num,m)=wo(num,m)+NNALPHA*hi(m)*d;
       end
       wo(num,3)=wo(num,3)+NNALPHA*(-1)*d;
       for N=1:HIDDENNO
           dj(N)=hi(N)*(1-hi(N))*wo(num,N)*(e(t,INPUTNO+num)-o(num,t))*o(num,t)*(1-o(num,t));
           for M=1:INPUTNO
               wh(N,M)=wh(N,M)+NNALPHA*e(t,M)*dj(N);
           end
           wh(N,3)=wh(N,3)+NNALPHA*(-1)*dj(N);
       end
            
    s=snext;
    end
    if s<cost
        cost=s;
    end
    ccc=(C-cost)^2;
    if ccc<0.001
        break
    end
    ans=qvalue;
end
      tooc=[5.7886 6.5569 7.2863 7.7683 8.2996 8.9049];      


