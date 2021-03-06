% Created on WED Mar 8th 13:25:45 2017
% Propose: Script1: plot spectrum data
%
% Enviroment: Matlab 2015b
% @auththor: kevin


%################################Part 1: Plot spectrum data############################################
%% 
% FM 88-108   TV 470�C798   GSM1800 UL(cell-phone to station) 1710�C1740
% CDMA2000 1920�C1935  

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
  clear;clc;close all;
 load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1710_1740.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=1710:0.025:1740;
data=dataLevel;
clear dataLevel;
display('GSM 1800 data loaded!!!');
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
  
 %################################Part 2: Energy detection############################################
% Not right 2017-3-17
 %  %% with simulation data
%  clear; clc; close all;
%  len=10000;
%  p=100;
%  p_fa=zeros(1,p); p_md=zeros(1,p);
%  for i=1:p
%      source_info{1}.mean=80;  source_info{1}.std=10;  source_info{1}.visit_time=25;
%      source_info{2}.mean=40;  source_info{2}.std=8;  source_info{2}.visit_time=38;
%      source_info{3}.mean=15;  source_info{3}.std=3;  source_info{3}.visit_time=40;
%      source_info{4}.mean=10;  source_info{4}.std=5;  source_info{4}.visit_time=25;
%      noise_info.mean=5; noise_info.std=2; 
%      [y,rcs]=Generate_simulation_dataset(source_info,noise_info,len);
%      [cut_point]=Recursive_oneside_hypthesis_testing(y, 100);
%      thres=cut_point(end);
%       ycs=y>thres;
%       err=ycs-rcs;
%       
%       err_fa=(err==1);
%       p_fa(i)=sum(err_fa)/length(rcs);
%       
%       err_md=(err==-1);
%       p_md(i)=sum(err_md)/length(rcs);
%       
%       err_all=sum(abs(ycs-rcs))/length(rcs);
%  end
%  sum(p_fa)/length(p_fa)
%  sum(p_md)/length(p_md)
%  figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
%  subplot(2,1,2); plot(p_md); xlabel('Times'); ylabel('Probability');title('Miss detection probability');
%  display('finished!!');
%  print('FA_MD','-dpng');

%% Generate simulation data and test with energy detection methods
clear; clc; close all;
len=1e3;
exp_time=1e2;
p_fa=zeros(1,exp_time); p_md=zeros(1,exp_time); %detection results of each time
rr=zeros(1,exp_time); %signal-noise ratio
for ii=1:exp_time
    [s,rcs,rr(1,ii)]=Generate_simulation_dataset_v2(60,10,len); %import noise parameters with total length
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

 sum(rr)/exp_time
 sum(p_fa)/length(p_fa)
 sum(p_d)/length(p_d)
 figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
 subplot(2,1,2); plot(p_d); xlabel('Times'); ylabel('Probability');title('Detection probability');


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
  %thresholding thresholding
 freqwidth=0.25; %1MHz for FM signal
 [s_cs,s_thres]=SlidingThresholding(data,freqix,freqwidth);
 figure(5);  subplot(4,1,1);  plot(s_thres); title('Adaptive thresholding');%plot(1:length(cut_point),cut_point); title('Thresholding spectrum');
 xlabel('Times'); ylabel('Threshold value');
 subplot(4,1,2); imagesc(pro_data); xlabel('Frequencies'); ylabel('Time slots'); title('Orignal spectrum sample');
 %direct threholding
 cs=pro_data>thres;
 subplot(4,1,3); imagesc(cs); xlabel('Frequencies'); ylabel('Time slots'); title('ROHT');
 subplot(4,1,4); imagesc(s_cs); xlabel('Frequencies'); ylabel('Time slots'); title('Sliding-ROHT');
 print('Figs/Adoptive thresholding','-dpng');
 a_cs=cs; %select which kind of cs data
 save('Datasets/ChannelStatus_88_108.mat','a_cs');
 tot_time=start_ti:step_ti:stop_ti;
 tot_time=datestr(tot_time,'yyyy-mm-dd HH:MM:SS')
 save('Datasets/Timeindex_88_108.mat','tot_time');
 display('channel analysis finished!!!');
 %% Calculate major index
 %SCR ҵ��ۼ���
 dc=sum(cs,2)/length(freqix); %average each row
 figure(6); subplot(4,1,1); plot(dc); title('duty cycle');
 save('ChannelStauseDataset/Duty_cycle_1710_1740.mat','dc');
 %% plot dc
 load('ChannelStauseDataset/Duty_cycle_1710_1740.mat');
 timeix=datetime(tot_time,'InputFormat','yyyy-MM-dd HH:mm:SS');
 figure(7); plot(timeix,dc); title=('Dutty cycle'); xlabel('Date'); ylabel('Ratio')
 print('Figs/DC','-dpng');
 
 %%  CVD
 cvd=zeros(size(freqix));
 for freq=1:length(freqix)
    duration=FindZerosBlock(cs(:,freq)');
    cvd(freq)=sum(duration)/length(duration);
 end
 %% mobile service utilization, MSU
 energysum=sum(pro_data,2);
 msu=diff(energysum,1); %first-order difference
 histogram(msu,'Normalization','probability');
 %################################Part 2: Channel status model ############################################
 %% load data
 clear;clc;close all;
 cs_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\ChannelStatus_1710_1740.mat';
 timeix_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\Timeindex_1710_1740.mat';
 load(cs_path);
 load(timeix_path);
 timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
 seq_mat=a_cs; seq_ix=dateStamp; clear a_cs; clear dateStamp; clear cs_path; clear timeix_path;
 %% get a sequence 
 yesfigure=0;
 startix=find(timeix==datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 stopix=find(timeix==datetime('2015-12-17 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 if yesfigure
    figure(1); 
    imagesc(seq_mat(startix:stopix, :)); title('GSM1800 UL in 2015-12-16');
 end
 seq=seq_mat(startix:stopix,120);
 %shape data to string
 seq=int16(seq(:,1)');
 %assemble a string
 seq_str='';
for i=1:length(seq)
    seq_str=strcat(seq_str,int2str(seq(i)));
end
 %% split original data into train-set and test-set
% seq_str='11111111111111111001111111111111111111111111111111';
alp = alphabet(seq_str);
split_point=length(seq_str)-10;
seq_train=seq_str(1:split_point);
seq_test=seq_str((split_point+1):end); 
%algorithm
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'PST'};;
params.ab_size = size(alp);
params.d = 10;  %maxmium order of VMM order %D
params.m = 2; %only for LZ-MS
params.s = 8; %only for LZ-MS
params.pMin = 0.006; %only for PST
params.alpha= 0; %only for PST
params.gamma = 0.0006; %only for PST
params.r = 1.05; %only for PST
params.vmmOrder = params.d; 
 
disp('---------------------------------------------------');
disp('working with AB={0,1}');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    %create a VMM through training 
    jVmm = vmm_create(map(alp, seq_train),  ALGS{i}, params);
    disp(sprintf('Pr(0 | 00) = %f', vmm_getPr(jVmm, map(alp,'0'), map(alp,'00'))));
    disp(sprintf('Pr(1 | 00) = %f', vmm_getPr(jVmm, map(alp,'1'), map(alp,'00'))));
    % calculates the length in bits of the  "compressed" representation of
    % seq.  -log[ Pr ( seq | jVmm) ]
    disp(sprintf('-lg(Pr(tar))=%f', vmm_logEval(jVmm,map(alp, seq_test))));
    disp('--------')
    disp(' ');
end

%################################Part 3: Spatial analysis############################################
%% Spatial
clear; clc; close all;

orig_path='D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\MultiSensorsDataset\\%s';
orig_filename='Dataset_%s_106800_106800.mat';

load(sprintf(orig_path,'SensorIds.mat'));
device_n=length(SensorIds); 
%path{1}=char('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\MultiSensorsDataset\\Dataset_2_101700_101700.mat');
%path{2}=char('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\MultiSensorsDataset\\Dataset_3_101700_101700.mat');

inital_flag=1; time_sp=0; data_sp=0;
locateion_sp=zeros(2,device_n);

for n=1:device_n
    %double string-number transfer to avoid name conflicts sucha as '02' and '2'
    t_filename=sprintf(orig_filename, num2str(str2num(SensorIds{n})));
    t_path=sprintf(orig_path, t_filename);
    load(t_path);
    %load time-stamp and inital data size in the first load
    if inital_flag==1
        [time_n,freq_n]=size(dataLevel);
        time_sp=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
        data_sp=zeros(time_n,device_n);
        inital_flag=0;
    end
    %new load data do not match
    [t_time,t_freq]=size(dataLevel);
    if t_time~=time_n || t_freq~=freq_n
        error('device data set do not match, exit!!!!');
    end
    data_sp(:,n)=dataLevel(:,1); %make sure the two vector matches
    locateion_sp(1,n)=deviceInfo.Longitude;
    locateion_sp(2,n)=deviceInfo.Latitude;
end

display('finish spatial data load!!!');
%% re-save loaded data
save_path='D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/%s';
save(sprintf(save_path, 'Spdata_1068_1068.mat'),'data_sp');
save(sprintf(save_path,'Sptime_1068_1068.mat'),'dateStamp');
save(sprintf(save_path,'Splocation_1068_1068.mat'),'locateion_sp');

%% pure spatial analysis
clear; clc; close all;
orign_path='D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/%s';
load(sprintf(orign_path, 'Spdata_1068_1068.mat'));
load(sprintf(orign_path, 'Spdata_1068_1068.mat'));
load(sprintf(orign_path,'Splocation_1068_1068.mat'));
%% 
data_sum=sum(data_sp,1,'omitnan')./1e6;
lon=locateion_sp(1,:); lat=locateion_sp(2,:);
plot(lon,lat,data_sum)



