% Created on WED Mar 17th 13:25:45 2017
% Propose: Plot scripts 
%
% Enviroment: Matlab 2015b
% @auththor: kevin

%% Plot sumarized spectrum figure

clear; clc; close all;
yes_figure=1;
%FM
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_88_108.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);
freqix=88:0.025:108;
data=dataLevel;
clear dateStamp; clear dataLevel;
display('Data loading finished!!!');

%%  GSM1800
  clear;clc;close all;
 load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1710_1740.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);
freqix=1710:0.025:1740;
%cut off
f1=(1730-1710)/0.025;
freqix=freqix(f1:end);
data=dataLevel(:,f1:end);
clear dataLevel;
display('GSM 1800 data loaded!!!');

%% TV
clear;clc;close all;
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_600_798.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=600:0.025:798;
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);

f1=(720-600)/0.025; f2=(750-600)/0.025;
freqix=freqix(:,f1:f2);
data=dataLevel(:,f1:f2);
clear dateStamp; clear dataLevel;
display('TV data loaded!!!');


 %% process missing data
 start_ti=timeix(1);
 step_ti=timeix(2)-timeix(1);
 stop_ti=timeix(end);
 
 ti_slots=(stop_ti-start_ti)/step_ti+1;
 pro_data=0;
 if ti_slots==length(data(:,1))
     display('no missed data');
     pro_data=data;
 %unequal to the therotical value, must miss some data
 else
     display('need to process missed data!!!');
     pro_data=zeros(ti_slots,length(freqix)); %prepare to process the missing data
     ti=start_ti; ixt=1; %index for time slots
     ix=1; %index for the row of dataset
     while ti<= stop_ti
         if  ti~=timeix(ix) %need special treatment
            pro_data(ixt,:)=NaN;
            ixt=ixt+1;
         else %just move the exist data to another matrix
            pro_data(ixt,:)=data(ix,:); 
            ixt=ixt+1; ix=ix+1;
         end
         ti=ti+step_ti;
     end
 end

 %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h_1=figure(1); set(gca,'Position',[.05 .05 .9 .9]);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (FM 88-108MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HHPM')
  print('Figs/fm','-dpng','-r500');
  display('FM analysised...');
  
  %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h2=figure(3);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (GSM1800 1730-1740MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  axis([1730 1740 0 1.0]);
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HHPM')
  print('Figs/gsm','-dpng','-r500');
  display('GSM analysised...');
  
  %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h2=figure(3);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (TV 600-800MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  axis([720 750 0 1.0]);
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HH')
  print('Figs/tv','-dpng','-r500');
  display('TV analysised...');
  
  %% ROHT-ED
  n_data(isnan(n_data(:,1)),:)=[]; %elimate NaN elements
  [cut_point]=Recursive_oneside_hypthesis_testing(n_data, 100);
  thres=cut_point(end);
  %thresholding thresholding
 freqwidth=0.25; %1MHz for FM signal
 [s_cs,s_thres]=SlidingThresholding(n_data,freqix,freqwidth);
 figure(3);  
 %direct threholding
 cs=n_data>thres;
 s_cs=cs|s_cs;
 subplot(3,1,1); imagesc(freqix,time,n_data); xlabel('频率/MHz'); ylabel('时间'); title('频谱监测数据'); datetick('y','HHPM');
 subplot(3,1,2); imagesc(freqix,time,cs); xlabel('频率/MHz'); ylabel('时间'); title('ROHT-能量检测结果'); datetick('y','HHPM');
%  subplot(4,1,3);  plot(s_thres); %title('Adaptive thresholding');%plot(1:length(cut_point),cut_point); title('Thresholding spectrum');
%  title('Sliding-ROHT的阈值');xlabel('窗口（窗宽250kHz）'); ylabel('阈值（已归一化）'); axis([0 81 0 1.5]);
 subplot(3,1,3); imagesc(freqix,time,s_cs); xlabel('频率/MHz'); ylabel('时间'); title('Sliding-ROHT-能量检测结果（窗宽250kHz）'); datetick('y','HHPM');
 %print('Figs/Adoptive thresholding','-dpng');
  
 %% 
 data=Normalize_matrix(pro_data,0,1)>thres;
 time=datenum(timeix(1):step_ti:timeix(end));
 subplot(3,1,1); imagesc(freqix,time,data); xlabel('频率/MHz'); ylabel('时间'); title('频谱占用状态'); datetick('y','HHPM');
 scr=sum(data,2)/length(freqix); %average each row
 figure(1); subplot(3,1,2); plot(time(1:1000),scr(1:1000)); title('业务聚集度'); datetick('x',' mm/dd HH','keepticks'); xlabel('月/天 时'); ylabel('聚集度值');
 
 subplot(3,1,3);
 cvd=zeros(size(freqix));
 for freq=1:length(freqix)
    duration=FindZerosBlock(data(:,freq)');
    cvd(freq)=sum(duration)/length(duration);
 end
 cvd_ix=find(cvd>300);
 cvd(cvd_ix)=[];
 histogram(cvd,200,'Normalization','probability'); xlabel('空闲窗口数'); ylabel('频率');
 
 
  %% Plot simulation data 
t = 0;%time for spectrum access event happen
n = 0;%access times
a = [];%access- time table 
len=1e3;
traffic_data =zeros(1,len);
while(t<len)
    e = exprnd(80);%interval of spectrum access--following expontiental distribution
    t = t + e;
    n = n + 1;
    a = [a;n,t];
end    
a = a'; 

occupy = geornd(0.2,1,n);%following geo-distribution
%set average occupancy time 
for j = 1:n
    i = ceil(a(2,j));
    m = occupy(j);
    traffic_data(i:i+m) = 1;
end 


%make the start point number equal with  stop point
if traffic_data(end)>eps 
    traffic_data=[traffic_data 0.0];
end
traffic=traffic_data>eps;
start_point=find(diff(traffic)==1)+1;
stop_point=find((diff(traffic)==-1));

%randomly pick signal distributions
s_mean_mu=40; s_mean_theata=10;
s_var_mu=10; s_var_theata=3;
s_num=4;
noise_mu=10; noise_theate=3;

mu=random('Normal', s_mean_mu, s_mean_theata, 1, s_num);
theata=random('Normal', s_var_mu, s_var_theata, 1, s_num);
s_ix=floor(s_num*rand(length(start_point),1)+1);
ns=random('Normal',noise_mu,noise_theate,1,length(traffic_data));


figure(1); set(gca,'Position',[.05 .05 .9 .9]); subplot(2,1,1);stem(traffic_data); xlabel('时间'); ylabel('占用状态'); title('（a）频谱接入状态');
axis([0 len eps 1.2]);
hold on;
for i=1:length(s_ix)
    text(start_point(i),1.1,sprintf('PU-%s',int2str(s_ix(i))));
end

traffic_data=ns;

for i=1:length(start_point)
     start_ix=start_point(i);  stop_ix=stop_point(i);
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)+random...
         ('Normal',mu(s_ix(i)),theata(s_ix(i)),1,stop_ix-start_ix+1);  %access from a randomly picked first user
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)-ns(start_ix:stop_ix); %remove the noise part
end

subplot(2,1,2);plot(traffic_data); xlabel('时间'); ylabel('频谱能量'); title('（b）频谱能量变化');
axis([0 len -2*noise_mu 2*s_mean_mu]);

print('Figs/sim','-dpng','-r500');
%% ROHT Process
rcs=traffic;
s=traffic_data;
split=floor(length(s)/2);
train_s=s(1:split);
[cut_point]=Recursive_oneside_hypthesis_testing(train_s, 100);
test_s=s(split+1:end);
test_rcs=rcs(split+1:end);
thres=cut_point(end);
test_ycs=test_s>thres;
err=test_ycs-test_rcs;
      
err_fa=(err==1);
p_fa=sum(err_fa)/length(test_rcs);
      
err_d=(err==0);
p_d=sum(err_d)/length(test_rcs);



%% ROHT-ED with sim
clear; clc; close all;
len=1e3;
exp_time=1e3;
p_fa=zeros(1,exp_time); p_md=zeros(1,exp_time); %detection results of each time
rr=zeros(1,exp_time); %signal-noise ratio
for ii=1:exp_time
    [s,rcs,rr(1,ii)]=Generate_simulation_dataset_v2(42,3,len); %import noise parameters with total length
    split=floor(length(s)/2);
    train_s=s(1:split);
    [cut_point]=Recursive_oneside_hypthesis_testing(train_s, 100);
    
    test_s=s(split+1:end);
    test_rcs=rcs(split+1:end);
    thres=cut_point(end);
    test_ycs=test_s>thres;
    err=test_ycs-test_rcs;
      
    err_fa=(err==1);
    p_fa(ii)=sum(err_fa)/length(test_rcs);
      
    err_d=(err==0);
    p_d(ii)=sum(err_d)/length(test_rcs);
      
    if p_d(ii)<=0.6
          display('abnormal!!!');
    end
      
      err_all=sum(abs(test_ycs-test_rcs))/length(test_rcs);
end
 
 %sum(rr)/exp_time
 display('False Alarm:');
 sum(p_fa)/length(p_fa)
 min(p_fa),max(p_fa)
 display('Detected:');
 sum(p_d)/length(p_d)
 max(p_d),min(p_d)
 figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
 subplot(2,1,2); plot(p_d); xlabel('Times'); ylabel('Probability');title('Detection probability');

%% Plot energy-detection ROC curve
clear; clc; close all;
snr=[-8,-9,-10,-11];
num=length(snr);
figure(1);
ix=['（a）','（b）','（c）','（d）'];
color={'--g','--b','--k','--r'};
color2={'g','b','k','r'};

figure(1)
for i=1:length(snr)
    [Pf,Pd,Pd_the]=Energy_detection_simulation(snr(i));
    subplot(2,2,i);
    plot(Pf, Pd,color{i},'LineWidth',2,'MarkerEdgeColor',color2{i}, 'MarkerSize',6); hold on; %paractical
    plot(Pf, Pd_the, color2{i},'LineWidth',2); %theory
    title(sprintf('SNR=%s dB',snr(i)));
    xlabel({['P_{fa}'],[ix(i)]}); ylabel('P_{d}');
    legend('实验值','理论值');
end
xlabel('P_{fa}'); ylabel('P_{d}');
%legend('SNR=-8dB','SNR=-9dB','SNR=-10dB','SNR=-11dB');
axis([0 0.5 0.5 1]);

print('Figs/ed','-dpng','-r500');

%% 

  
  
