%% load data 
clear; clc; close all;
orign_path='D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/%s';
load(sprintf(orign_path, 'Spdata_1730_1740.mat'));
load(sprintf(orign_path, 'Sptime_1730_1740.mat'));
load(sprintf(orign_path,'Splocation_1730_1740.mat'));

%% load another data
clear; clc; close all;
orign_path='D:/Code/WorkSpace/ThesisCode/Src/5_Results/Datasets/%s';
load(sprintf(orign_path, 'Data.mat'));
load(sprintf(orign_path, 'Position(2015).mat'));
load(sprintf(orign_path, 'StatData.mat'));
value=zeros(1,60);
for i=1:60
    value(1,i)=StatData{i}(1);
end

%% 
%clear; clc; close all;



%%
% shape data
% x=locateion_sp(1,:)';
% y=locateion_sp(2,:)';
% x=Clip{1,1};
% y=Clip{1,2};
% z=Clip{1,3};
x=Lon; y=Lat; z=value';
points=[x,y]; len=60;
x=(x-min(x))/((max(x)-min(x)));
y=(y-min(y))/((max(y)-min(y)));
dis_nor=zeros(1,len*len);
dis_km_Haversine=zeros(size(dis_nor));
dis_km_Pythagoras=zeros(size(dis_nor));
for i=1:len
    for j=i:len
    dis_nor(1,i*j)=sqrt((points(i,1)-points(j,1))^2+(points(i,2)-points(j,2))^2);
    [dis_km_Haversine(1,i*j),dis_km_Pythagoras(1,i*j)]=lldistkm(points(i,:),points(j,:));
    end
end
% abnormal_ix=[];
% abnormal_ix=find(z<70);
% x(abnormal_ix)=[]; y(abnormal_ix)=[]; z(abnormal_ix)=[];

% z=(z-min(z))/((max(z)-min(z)));
% z=data_sp(255,:)';
[X,Y]=meshgrid(-0.01:0.005:1.01);
max_dis=1.4;

%% empircal variogram
 v = variogram([x y],z,'plotit',false,'maxdist',max_dis);
 
 %% fit theotical model
 figure(1);
 [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','gaussian','solver', 'fminsearchbnd');
 xlabel('归一化距离','FontSize',12); ylabel('变差值/dB\muVm^{-1}','FontSize',12);
 [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
%  print('Figs/varigram','-dpng','-r500');

figure(4);
subplot(2,2,1); [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','blinear','solver', 'fminsearchbnd');
 [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
imagesc(X(1,:),Y(:,1),sqrt(Zvar)); axis image; axis xy;
 h=colorbar; xlabel(h,'能量估计标准差/dB\muVm^{-1}','FontSize',12);
xlabel('相对经度','FontSize',12); ylabel('相对纬度','FontSize',12);

subplot(2,2,2); [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','circular','solver', 'fminsearchbnd');
 [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
imagesc(X(1,:),Y(:,1),sqrt(Zvar)); axis image; axis xy; 

subplot(2,2,3); [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','exponential','solver', 'fminsearchbnd');
 [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
imagesc(X(1,:),Y(:,1),sqrt(Zvar)); axis image; axis xy; 

subplot(2,2,4); [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable','solver', 'fminsearchbnd');
 [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
imagesc(X(1,:),Y(:,1),sqrt(Zvar)); axis xy; axis image;
%
path='D:/doc/PapaerLibrary/Figures/kriging_err';
print(path,'-dpng','-r500');
%% plot kriging
 figure(2);
 imagesc(X(1,:),Y(:,1),Zhat); axis image; axis xy; xlabel('相对经度','FontSize',12); ylabel('相对纬度','FontSize',12);
 h = colorbar; hold on; xlabel(h,'频谱能量/dB\muVm^{-1}','FontSize',12)
% print('Figs/kriging_result','-dpng','-r500');

figure(3);
imagesc(X(1,:),Y(:,1),sqrt(Zvar)); axis image; axis xy; xlabel('相对经度','FontSize',12); ylabel('相对纬度','FontSize',12);
h=colorbar; xlabel(h,'能量估计标准差/dB\muVm^{-1}','FontSize',12)
% print('Figs/kriging_error','-dpng','-r500');

%% ITU-P.R 1546
figure(4);
d=1:0.5:100;
f=1800; 
h1=100;
h2=5;
t=50;
tca=10; eff2=15;
for x=1:length(d)
    y(x)=P1546FieldStr(d(x),f,t,h1,h1)...
    -Step_12(f,tca)-Step_14(h2,f,'Land','open',d,h1);
end
yy=(y(1)-y)*(y(1)-y);
plot(d,yy,'LineWidth',1.0);
hold on;
[dum,dum,dum,vstruct] = variogramfit(v.distance,sqrt(v.val),[],[],[],'model','exponential','solver', 'fminsearchbnd');
dis=v.dis;
plot()
xlabel('距离/km','FontSize',12); ylabel('频谱能量衰减值/dB\muVm^{-1}','FontSize',12);
%print('Figs/itu.p.r.1546','-dpng','-r500');
%% Variogram with progration
points=[x,y];
len=60;
d1=zeros(1,factorial(len)); d2=zeros(1,factorial(n));
for i=1:60
    for j=i:60
        [d1(1,i*j),d2(1,i*j)]=lldistkm;
end
