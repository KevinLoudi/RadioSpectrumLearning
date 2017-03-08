% Created on WED Mar 8th 13:25:45 2017
% Propose: Script1: plot spectrum data
%
% Enviroment: Matlab 2015b
% @auththor: kevin


%################################Part 1: Plot spectrum data############################################
%% 
% FM 88-108   TV 470¨C798   GSM1800 UL(cell-phone to station) 1710¨C1740
% CDMA2000 1920¨C1935  

clear; clc; close all;
%FM
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_88_108.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=88:0.025:108;
data=dataLevel;
clear dateStamp; clear dataLevel;

%% 
figure(1);
  subplot(2,1,1); plot(freqix,data(1,:)); title('FM Spectrum in frequency');
  subplot(2,1,2); imagesc(data); title('FM Spectrum in time and frequency'); 
 print('Figs/FM','-dpng');
 display('FM analysised...');
  
  %% 
 %TV
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_600_798.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=600:0.025:798;
data=dataLevel;
clear dateStamp; clear dataLevel;

%% 
figure(2);
  subplot(2,1,1); plot(freqix,data(1,:)); title('TV Spectrum in frequency');
  subplot(2,1,2); imagesc(data); title('TV Spectrum in time and frequency'); 
  print('Figs/TV','-dpng');
  display('TV analysised...');
  
  %%  GSM1800
 load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1710_1740.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=1710:0.025:1740;
data=dataLevel;
clear dateStamp; clear dataLevel;
%% 
figure(3);
  subplot(2,1,1); plot(freqix,data(1,:)); title('GSM1800 UL Spectrum in frequency');
  subplot(2,1,2); imagesc(data); title('GSM1800 UL Spectrum in time and frequency'); 
  print('Figs/GSM1800','-dpng');
 display('GSM1800 analysised...');
  
  %% CDMA2000
   load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1920_1935.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=1920:0.025:1935;
data=dataLevel;
clear dateStamp; clear dataLevel;
%% 
figure(4);
  subplot(2,1,1); plot(freqix,data(1,:)); title(' CDMA2000 UL Spectrum in frequency');
  subplot(2,1,2); imagesc(data); title(' CDMA2000 UL Spectrum in time and frequency'); 
  print('Figs/CDMA2000','-dpng');
  display('CDMA2000 analysised...');
  display('All data analysising finished!!!');
  
 %################################Part 1: Energy detection############################################
 %% with simulation data
 clear; clc; close all;
 len=10000;
 p=1000;
 p_fa=zeros(1,p); p_md=zeros(1,p);
 for i=1:p
     source_info{1}.mean=50;  source_info{1}.std=10;  source_info{1}.visit_time=3;
     source_info{2}.mean=30;  source_info{2}.std=8;  source_info{2}.visit_time=8;
     source_info{3}.mean=15;  source_info{3}.std=3;  source_info{3}.visit_time=10;
     source_info{4}.mean=10;  source_info{4}.std=5;  source_info{4}.visit_time=6;
     noise_info.mean=5; noise_info.std=2; 
     [y,rcs]=Generate_simulation_dataset(source_info,noise_info,len);
     [cut_point]=Recursive_oneside_hypthesis_testing(y, 100);
     thres=cut_point(end);
      ycs=y>thres;
      err=ycs-rcs;
      
      err_fa=(err==1);
      p_fa(i)=sum(err_fa)/length(rcs);
      
      err_md=(err==-1);
      p_md(i)=sum(err_md)/length(rcs);
      
      err_all=sum(abs(ycs-rcs))/length(rcs);
 end
 sum(p_fa)/length(p_fa)
 sum(p_md)/length(p_md)
 figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
  subplot(2,1,2); plot(p_md); xlabel('Times'); ylabel('Probability');title('Miss detection probability');
 display('finished!!');
 print('FA_MD','-dpng');
 %% thresholding
 [cut_point]=Recursive_oneside_hypthesis_testing(y, 100);
 thres=cut_point(end);
 %% make decision on channel status
 ycs=y>thres;
 err=sum(abs(ycs-rcs))/length(rcs);

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
 
 %% energy detction with real data
 [cut_point]=Recursive_oneside_hypthesis_testing(data, 100);
 thres=cut_point(end);
 figure(5); subplot(4,1,1);  plot(s_thres); title('Adaptive thresholding');%plot(1:length(cut_point),cut_point); title('Thresholding spectrum');
 xlabel('Times'); ylabel('Threshold value');
 subplot(4,1,2); imagesc(pro_data);
 %direct threholding
 cs=pro_data>thres;
 subplot(4,1,3); imagesc(cs);
 %thresholding thresholding
 freqwidth=3; %1MHz for FM signal
 [s_cs,s_thres]=SlidingThresholding(data,freqix,freqwidth);
 subplot(4,1,4); imagesc(s_cs);
 
 %% Calculate major index
 %dc
 dc=sum(cs,2)/length(freqix); %average each row
 figure(6); subplot(4,1,1); plot(start_ti:step_ti:stop_ti,dc); title('duty cycle');
 %%  cvd
 cvd=zeros(size(freqix));
 for freq=1:length(freqix)
    duration=FindZerosBlock(cs(:,freq)');
    cvd(freq)=sum(duration)/length(duration);
 end
 

 

 
 
 
 
 