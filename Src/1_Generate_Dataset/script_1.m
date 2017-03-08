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
 time_now=datetime(now);
 
 %% with FM data
 start_ti=timeix(1)
 step_ti=timeix(2)-timeix(1)
 stop_ti=timeix(end)
 
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
 
 
%  while ti<=stop_ti
%      
%  end
 
 