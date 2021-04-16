clc
clear all
%{
MEC 服务器的 CPU 工作频率 4GHz
任务的传输数据大小 在 300~1200KB 之间随机分布
任务所需的 CPU 工作周期数 0.1*10^9~1*10^9 之间随机分布
任务的最大完成时间 在 0.5~5s 之间随机分布
移动设备的本地 CPU 工作频率 0.2~1GHz
信道带宽 0.2MHz
 x 应该是随机产生的 x（1，：）是数据长度  x（2，：）是需要的工作周期数
                   x(3,:)是最大延迟  x(4,:)是信道参数
%}
x=[735534.5  652180.2  800773.8  715923.5  921134.9  892048.3  743773.8  790930.5  705540.3  869165.0  715923.5  921134.9  892048.3  743773.8  787157.0 984575.0 702192.0 768663.0 657394.0
   508229664 430814447 608749645 600559994 561337651 453478757 598871791 466361361 527207144 481861654 600559994 561337651 453478757 598871791 461197222 357248933 491403999 397368822 786297177
   3.052794  2.430809  2.423213  2.735298  2.517545  3.368810  2.838018  2.814051  2.229461  2.322614  2.735298  2.517545  3.368810  2.838018  1.802526 2.704607 3.348886 3.30791 2.217739
   2.925414  6.147478  3.530716  1.833023  2.653565  7.568968  5.010552  35.71505  4.165090  2.917472  1.833023  2.653565  7.568968  5.010552  3.925414 , 6.147478 , 3.530716 , 2.833023 , 3.653565];
%n(:,1)=(wgn(1,1000000,0));%信道噪声
%hy=[n n n n n n n n n n n n n n n n n n n n ];%按理来说每个信道的噪声是不同的
%y=(hy);
t=5;
k=10^(-26);
for nn=1:15;%loop for different task number
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
%MEC time and energy（no noise）
for n=1:nn
       m=1;
        for pi=0.001:0.0005:10^-0.7
            ri=0.2*10^6*log2(1+pi*x(4,n)^2);%upload transform speed
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
            ri=0.2*10^6*log2(1+pi*x(4,n)^2);
            ec(n,m)=pi*x(1,n)/ri;
            m=m+1;
        end
    
end
GC=0.2*tc+0.8*ec;
for n=1:nn
gc(n)=min(GC(n,:));
end
g=[gl;gc];
for n=1:nn
G(n)=min(g(:,n));
end
C(nn)=sum(G) %total cost （no noise）

%MEC time and energy（noise）
n=1;
for l=1:1000
for n=1:nn
       m=1;
        for pi=0.001:0.00099:0.2
            ri=0.2*10^6*log2(1+pi*abs(x(4,n)+y(l,n))^2);
            tcn(n,m)=x(1,n)/ri+x(2,n)/(4*10^9);
           if tcn(n,m)>x(3,n)
                tcn(n,m)=100;
           end
            m=m+1;
            
        end
    
end
for n=1:nn
       m=1;
        for pi=0.001:0.00099:0.2
            ri=0.2*10^6*log2(1+pi*abs(x(4,n)+y(l,n))^2);
            ecn(n,m)=pi*x(1,n)/ri;
            m=m+1;
        end
    
end

GCn=0.2*tcn+0.8*ecn;
for n=1:nn
gcn(l,n)=min(GCn(n,:));
end
end
for n=1:nn
gccn(n)=mean(gcn(:,n));%expectation
end

%MEC minimum total cost
for n=1:nn
[rowc(n) colc(n)]=find(GCn(n,:)==min(GCn(n,:)));
end
for n=1:nn
TCn(n)=tcn(n,colc(n));
end
for n=1:nn
ECn(n)=ecn(n,colc(n));
end

%total cost （noise）
g=[gl;gccn];
for n=1:nn
Gn(n)=min(g(:,n));
end
Cn(nn)=sum(Gn)
end
plot(Cn)
hold on
plot(C)