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
  [X,Y]=meshgrid(freqix,time);
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
  
