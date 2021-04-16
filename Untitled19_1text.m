clc
clear all
GENMAX=200000;
STATION=2;
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
x=[787157.0 984575.0 702192.0 768663.0 657394.0 715315.0 872890.0 780967.0 591001.0 729557.0 764296.0 804291.0 1057499.0 800046.0 782501.0 735534.5 652180.2 800773.8 715923.5 921134.9 892048.3 743773.8 790930.5 705540.3 869165.0
   461197222 357248933 491403999 397368822 786297177 597629733 519093522 589905788 628979355 516016799 576636122 430706177 631546533 561690777 598172211 508229664 430814447 608749645 600559994 561337651 453478757 598871791 466361361 527207144 481861654
   1.802526 2.704607 3.348886 3.30791 2.217739 2.039008 3.435385 3.056378 2.195885 2.287386    3.432247 3.727381 3.738657 3.433317 2.427195 3.052794 2.430809 2.423213 2.735298 2.517545 3.368810 2.838018 2.814051 2.229461 2.322614
   ];
hi=[ 3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565, 7.568968 , 5.010552 , 35.71505 , 4.165090 , 2.917472;
          ];
      
      I=zeros(1,10);

n(:,1)=(wgn(1,1000000,0));%信道噪声
hy=[n n n n n n n n n n];%按理来说每个信道的噪声是不同的
y=(hy);

for n=1:nn
       m=1;
        for fi=0.2*10^9:0.001*10^8:1*10^9
            t(n,m)=x(2,n)/fi;
            if t(n,m)>x(3,n)
                t(n,m)=100;
            end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for fi=0.2*10^9:0.001*10^8:1*10^9
            ei(n,m)=x(2,n)*fi^2*10^(-26);
            m=m+1;
        end
    
end
GL=0.2*t+0.8*ei;
for n=1:nn
gl(n)=min(GL(n,:));
end

for n=1:nn
    for s=1:10
       m=1;
        for pi=0.001:0.0005:10^-0.7
             for lr=1:100
            r(lr)=0.2*10^6*log2(1+pi*(hi(s)+y(lr,s))^2/(I(s)+1));
            end
            ri=mean(r);
            tc(n,m,s)=x(1,n)/ri+x(2,n)/(4*10^9);
           if tc(n,m,s)>x(3,n)
                tc(n,m,s)=100;
           end
            m=m+1;          
        end
    end

    for s=1:10
       m=1;
        for pi=0.001:0.0005:10^-0.7
            for lr=1:100
            r(lr)=0.2*10^6*log2(1+pi*(hi(s)+y(lr,s))^2/(I(s)+1));
            end
            ri=mean(r);
            ec(n,m,s)=pi*x(1,n)/ri;
            m=m+1;
        end
    end

GC=0.2*tc+0.8*ec;

    aa(n)=1000000;
    for s=1:10
gc(n)=min(GC(n,:,s));
if gc(n)<aa(n)
    aa(n)=gc(n);
    num=s;
[min_val, position_min] = min(GC(n,:,s)); 
[p,q] = ind2sub(size(GC(n,:,s)),position_min);
end
    end
g=[gl(n);aa(n)];
G(n)=min(g);
if g(1,1)<g(2,1)
    sa(n)=0;
else
    sa(n)=1;
    I(num)=(pi+q*0.0005)*hi(num)^2+I(num);%干扰更新
end
end
C=sum(G);
cl=sum(sum(gl))
cost=cl;
cc=sum(C)
for n=1:nn
if gl(n)<aa(n)
    e(n,1)=0;
    e(n,2)=1;
    e(n,3)=0;
else
    e(n,1)=1;
    e(n,2)=0;
    e(n,3)=1;
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
i=1;j=1;sa=0;
for loop=1:GENMAX
    s=0;
    m=1;
    I=zeros(1,10);
    for t=1:LEVEL
        ran=rand();
       
        if ran<EPSILON
        
        a= randi(2);
        crl=inf;
               for fi=0.2*10^9:0.001*10^8:1*10^9
                   til=x(2,t)/fi;
                   if til>x(3,t)
                      til=inf;
                    end
                   eil=x(2,t)*fi^2*10^(-26);
                   gil=0.2*til+0.8*eil;
                   
                   if gil<crl
                       crl=gil;
                   end
               end
                crc=inf;
               for sa=1:10;
               for pi=0.001:0.0005:10^-0.7
                 for l=1:100
                    rir(l)=0.2*10^6*log2(1+pi*(hi(sa)+y(l,sa))^2/(I(sa)+1));
                end
                    ric=mean(rir);
                 
               tic=x(1,t)/ric+x(2,t)/(4*10^9);
               if tic>x(3,t)
                      tic=inf;
               end
               eic=pi*x(1,t)/ric;
               gic=0.2*tic+0.8*eic;
          
                   if gic<crc
                       crc=aa(t);
                       pim=pi;
                       sac=sa;
                   end
               end
               end
           if a==1
               c=crl;
           else
               c=crc;
           end
           if crl< crc
               rc=crl;
           else
               rc= crc;
           end         
        else
            for i=1:HIDDENNO
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
            
            a=round(o(num,t)+1);
             cfl=inf;
             for fi=0.2*10^9:0.001*10^8:1*10^9
                   til=x(2,t)/fi;
                   if til>x(3,t)
                      til=inf;
                    end
                   eil=x(2,t)*fi^2*10^(-26);
                   gil=0.2*til+0.8*eil;
                  
                   if gil<cfl
                       cfl=gil;
                   end
             end
               cfc=inf;
               for sa=1:10;
                for pi=0.001:0.0005:10^-0.7
                 for l=1:100
                    rir(l)=0.2*10^6*log2(1+pi*(hi(sa)+y(l,sa))^2/(I(sa)+1));
                  end
                    ric=mean(rir);
                 
               tic=x(1,t)/ric+x(2,t)/(4*10^9);
               if tic>x(3,t)
                      tic=inf;
               end
               eic=pi*x(1,t)/ric;
               gic=0.2*tic+0.8*eic;
          
                   if gic<cfc
                       cfc=aa(t);
                       pim=pi;
                       sac=sa;
                   end
                end
               end
           if a==1
               c=cfl;
           else
               c=cfc;
           end
           if cfl< cfc
               rc=cfl;
           else
               rc=cfc;
           end    
                
        end
          if a==1
          I(sac)=I(sac);
          else
          I(sac)=pim*hi(sac)^2+I(sac);
          end
        snext=s+c;
       %[rownext colnext]=find(qvalue(t+1,:)==max(qvalue(t+1,:)));
        %anext=colnext;
        REWARD=(1-c/cl);
        if c==rc
        qv=qvalue(t,a)+ALPHA*(10*REWARD-qvalue(t,a));
        else
        qv=qvalue(t,a)+ALPHA*(REWARD-qvalue(t,a));
        end
        
    qvalue(t,a)=qv;
       for i=1:HIDDENNO
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
    ccc=(cc-cost)^2;
    if ccc<0.001
        break
    end
    ans=qvalue;
end
            

