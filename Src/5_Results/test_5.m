points=[Lon,Lat];
len=60;
d1=zeros(1,len*len); d2=zeros(1,len*len);
for i=1:60
    for j=i:60
        [d1(1,i*j),d2(1,i*j)]=lldistkm(points(i,:),points(j,:));
        var(1,i*j)=abs(value(i)-value(j));
    end
end

steps=25;
dmax=max(d1);
dstep=0:dmax/(steps-1):dmax;
vario=zeros(1,steps-1);

for i=1:(steps-1)
    %find samples within a certain range
    ix=find(d1<dstep(i+1));
    ix=find(d1(ix)>=dstep(i));
    %average variation
    vario(1,i)=sum(var(ix))/length(ix);
end

figure(4);
mark={'-','-.','--','-+'};
subplot(2,2,1);
h1_r=[10, 50, 100, 200];
for jj=1:4
 d=1:0.1:20;
 f=100; 
 h1=h1_r(jj);
 h2=5;
 t=50;
 tca=10; eff2=15;
 for x=1:length(d)
    y(x)=P1546FieldStr(d(x),f,t,h1,h1)...
    -Step_12(f,tca)-Step_14(h2,f,'Land','open',d,h1);
 end
 yy=y;%y(1)-y;
 plot(d,yy(1:length(d)),mark{jj},'LineWidth',1.0); hold on;

end
legend('发射天线高 h = 10m','发射天线高 h = 50m','发射天线高 h = 100m','发射天线高 h = 200m')

subplot(2,2,2);
h2_r=[1, 5, 10, 20];
for jj=1:4
 d=1:0.1:20;
 f=100; 
 h1=100;
 h2=h2_r(jj);
 t=50;
 tca=10; eff2=15;
 for x=1:length(d)
    y(x)=P1546FieldStr(d(x),f,t,h1,h1)...
    -Step_12(f,tca)-Step_14(h2,f,'Land','open',d,h1);
 end
 yy=y;%y(1)-y;
 plot(d,yy(1:length(d)),mark{jj},'LineWidth',1.0); hold on;

end
legend('接收天线高 h = 1m','接收天线高 h = 5m','接收天线高 h = 10m','接收天线高 h = 20m');

subplot(2,2,3);
f_r=[100, 900, 1700, 2600];
for jj=1:4
 d=1:0.1:20;
 f=f_r(jj); 
 h1=100;
 h2=5;
 t=50;
 tca=10; eff2=15;
 for x=1:length(d)
    y(x)=P1546FieldStr(d(x),f,t,h1,h1)...
    -Step_12(f,tca)-Step_14(h2,f,'Land','open',d,h1);
 end
 yy=y;%y(1)-y;
 plot(d,yy(1:length(d)),mark{jj},'LineWidth',1.0); hold on;

end
legend('发射频率 f = 100MHz','发射频率 f = 900MHz','发射频率 f = 1700MHz','发射频率 f = 2600MHz')

subplot(2,2,4);
field={'open', 'urban', 'dense', 'land'};
for jj=1:4
 d=1:0.1:20;
 f=100; 
 h1=100;
 h2=5;
 t=50;
 tca=10; eff2=15;
 for x=1:length(d)
    y(x)=P1546FieldStr(d(x),f,t,h1,h1)...
    -Step_12(f,tca)-Step_14(h2,f, 'Land',field{jj},d,h1);
 end
 yy=y;%y(1)-y;
 plot(d,yy(1:length(d)),mark{jj},'LineWidth',1.0); hold on;

end
legend('传播环境：空旷地','传播环境：城市','传播环境：大城市','传播环境：郊区')
xlabel('距离/km','FontSize',12); ylabel('频谱能量衰减/dB\muVm^{-1}','FontSize',12);
path='D:/doc/PapaerLibrary/Figures/itu_1546';
print(path,'-dpng','-r500');
