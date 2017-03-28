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
legend('�������߸� h = 10m','�������߸� h = 50m','�������߸� h = 100m','�������߸� h = 200m')

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
legend('�������߸� h = 1m','�������߸� h = 5m','�������߸� h = 10m','�������߸� h = 20m');

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
legend('����Ƶ�� f = 100MHz','����Ƶ�� f = 900MHz','����Ƶ�� f = 1700MHz','����Ƶ�� f = 2600MHz')

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
legend('�����������տ���','��������������','���������������','��������������')
xlabel('����/km','FontSize',12); ylabel('Ƶ������˥��/dB\muVm^{-1}','FontSize',12);
path='D:/doc/PapaerLibrary/Figures/itu_1546';
print(path,'-dpng','-r500');
