clc
clear all
GENMAX=2000000;
STATION=20000;
ACTIONNO=2000;
ALPHA=0.1;
GAMMA=0.9;
EPSILON=0.3;
  nn=25;
GOAL=55;
LEVEL=nn;
INPUTNO=20000;
HIDDENNO=2000;
OUTPUTNO=2000;
NNALPHA=3;
count=0;
qvalue=0.999*rand(nn,10,2000);
x(:,:,1)=[787157.0 984575.0 702192.0 768663.0 657394.0
          461197222 357248933 491403999 397368822 786297177
   1.802526 2.704607 3.348886 3.30791 2.217739];
x(:,:,2)=[715315.0 872890.0 780967.0 591001.0 729557.0
   597629733 519093522 589905788 628979355 516016799
  2.039008 3.435385 3.056378 2.195885 2.287386
   ];
x(:,:,3)=[764296.0 804291.0 1057499.0 800046.0 782501.0
   576636122 430706177 631546533 561690777 598172211
   3.432247 3.727381 3.738657 3.433317 2.427195
   ];
x(:,:,4)=[735534.5  652180.2  800773.8  715923.5  921134.9
          508229664 430814447 608749645 600559994 561337651
          3.052794  2.430809  2.423213  2.735298  2.517545
          ];
x(:,:,5)=[892048.3  743773.8  790930.5  705540.3  869165.0
          453478757 598871791 466361361 527207144 481861654
          3.368810  2.838018  2.814051  2.229461  2.322614
          ];

    
        hi=[ 3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565, 7.568968 , 5.010552 , 35.71505 , 4.165090 , 2.917472;
          ];
    
      I=zeros(1,10);
n(:,1)=(wgn(1,1000000,0));%信道噪声
hy=[n n n n n n n n n n];%按理来说每个信道的噪声是不同的
y=(hy);
for l=1:nn/5
for n=1:5
       m=1;
        for fi=0.2*10^9:(0.8*1e+9)/999:1*10^9
            t(n,m,l)=x(2,n,l)/fi;
            if t(n,m,l)>x(3,n,l)
                t(n,m,l)=100;
            end
            m=m+1;
            
        end
    
end
for n=1:5
       m=1;
        for fi=0.2*10^9:(0.8*1e+9)/999:1*10^9
            e(n,m,l)=x(2,n,l)*fi^2*10^(-26);
            m=m+1;
        end
    
end
GL=0.2*t+0.8*e;
for n=1:5
gl(l,n)=min(GL(n,:,l));
end
for n=1:5
[rowl(n,l) coll(n,l)]=find(GL(n,:,l)==min(GL(n,:,l)));

TL(n,l)=t(n,coll(n,l));

EL(n,l)=e(n,coll(n,l));
end
for n=1:5
    for s=1:10
       m=1;
        for pi=0.001:0.1985/999:10^-0.7
             for lr=1:100
            r(lr)=0.2*10^6*log2(1+pi*(hi(s))^2/(I(s)+1));
            end
            ri=mean(r);
            tc(n,m,s,l)=x(1,n,l)/ri+x(2,n,l)/(4*10^9);
           if tc(n,m,s,l)>x(3,n,l)
                tc(n,m,s,l)=100;
           end
            m=m+1;
            
        end
    end

    for s=1:10
       m=1;
        for pi=0.001:0.1985/999:10^-0.7
            for lr=1:100
            r(lr)=0.2*10^6*log2(1+pi*(hi(s))^2/(I(s)+1));
            end
            ri=mean(r);
            ec(n,m,s,l)=pi*x(1,n,l)/ri;
            m=m+1;
        end
    end

GC=0.2*tc+0.8*ec;

    aa(l,n)=1000000;
    for s=1:10
gc(l,n)=min(GC(n,:,s,l));
if gc(l,n)<aa(l,n)
    aa(l,n)=gc(l,n);
    num=s;
[min_val, position_min] = min(GC(n,:,s,l)); 
[p,q(n,l)] = ind2sub(size(GC(n,:,s,l)),position_min);
end
    end
g=[gl(l,n);aa(l,n)];
G(n,l)=min(g);
if g(1,1)<g(2,1)
    sa(l,n)=0;
    qvalue((l-1)*5+n,num,coll(n,l))=1;
else
    sa(l,n)=1;
    I(num)=(0.001+q(n,l)*0.1985/999)*hi(num)^2+I(num);%干扰更新
    qvalue((l-1)*5+n,num,q(n,l)+1000)=1;
end
end
C(l)=sum(G(:,l));
end
cl=sum(sum(gl))
cost=cl;
cc=sum(C)
for i=1:HIDDENNO
    for j=1:INPUTNO+1
        wh(i,j)=1-2*rand();
    end
end
i=1;j=1;
for i=1:OUTPUTNO
    for j=1:HIDDENNO
        wo(i,j)=1-2*rand();
    end
end
i=1;j=1;
for i=1:GENMAX
    sa=0;
    m=1;
    I=zeros(1,10);
    for t=1:LEVEL
     
        ran=rand();
        if ran<EPSILON
            s=randi(10);    
            a= randi(2000);
             if a<=1000
               fri=(0.8*rand()+0.2)*(1e+9);
               tri=x(2,t)/fri;
               eri=x(2,t)*fri^2*10^(-26);
               c=0.2*tri+0.8*eri;
             else
               pri=0.1995*rand();
               for l=1:100
               r(l)=0.2*10^6*log2(1+pri*(hi(s)+y(l,s))^2/(I(s)+1));
               end
               ri=mean(r);
               tci=x(1,t)/ri+x(2,t)/(4*10^9);
               eci=pri*x(1,t)/ri;
               c=0.2*tci+0.8*eci;
              
               end
                   
        else
            
             [max_val, position_max] = max(max(qvalue(1,:,:))); 
             [mobile rowl coll] = ind2sub(size(qvalue(1,:,:)),position_max);
            
            a=coll;
            s=rowl;
            if a<=1000
               fri=0.2*1e+9+coll*(0.8*1e+9)/999;
               tri=x(2,t)/fri;
               eri=x(2,t)*fri^2*10^(-26);
               c=0.2*tri+0.8*eri;
            else
                pri=0.001+(coll-1000)*0.1985/999;
                for l=1:100
                r(l)=0.2*10^6*log2(1+pri*(hi(s)+y(l,s))^2/(I(s)+1));
                end
                ri=mean(r);
               tci=x(1,t)/ri+x(2,t)/(4*10^9);
               eci=pri*x(1,t)/ri;
               c=0.2*tci+0.8*eci;
                
               end
        end
        snext=sa+c;
       %[rownext colnext]=find(qvalue(t+1,:)==max(qvalue(t+1,:)));
        %anext=colnext;
        
        REWARD=(1-c/cl);
        erro=(c-G(t))^2;
        if erro<=0.001
        qv=qvalue(t,s,a)+ALPHA*(3*REWARD-qvalue(t,s,a));
        else
        qv=qvalue(t,s,a)+ALPHA*(REWARD-qvalue(t,s,a));
        end
    qvalue(t,s,a)=qv;
    sa=snext;
    if a>1000
    I(s)=pri*hi(s)^2+I(s);
    else
    I(s)=I(s);
    end
    end
    if sa<cost
        cost=sa;
    end
    ccc=(cc-cost)^2;
    if ccc<0.01
        break
    end
    ans=qvalue;
    end
            

