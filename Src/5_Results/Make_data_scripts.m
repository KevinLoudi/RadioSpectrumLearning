%% Spatial
clear; clc; close all;

orig_path='D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\MultiSensorsDataset\\%s';
orig_filename='Dataset_%s_1730000_1740000.mat';

load(sprintf(orig_path,'SensorIds.mat'));
device_n=length(SensorIds); 

inital_flag=1; time_sp=0; 
locateion_sp=zeros(2,device_n);

%uniform data set range
 start_ti=datetime('2015-12-15 00:00:00','InputFormat','yyyy-MM-dd HH:mm:SS');
 step_ti=minutes(5); %5 minutes duration between 2 time slots
 stop_ti=datetime('2015-12-22 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS');
 freqix=1730000:25:1740000;
 exp_ti_slots=(stop_ti-start_ti)/step_ti+1; %expected time slots
 data_sp=zeros(exp_ti_slots,device_n);

for n=1:device_n
   
    %double string-number transfer to avoid name conflicts sucha as '02' and '2'
    t_filename=sprintf(orig_filename, num2str(str2num(SensorIds{n})));
    t_path=sprintf(orig_path, t_filename);
    load(t_path);
    timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
    %% process missing data
      pro_data=0;
     %deal with missing data
     if exp_ti_slots==length(dataLevel(:,1))
       display('no missed data');
       pro_data=dataLevel;
    %unequal to the therotical value, must miss some data
     else
    % display('need to process missed data!!!');
     pro_data=zeros(exp_ti_slots,length(freqix)); %prepare to process the missing data
     ti=start_ti; ixt=1; %index for time slots
     ix=1; %index for the row of dataset
     
     %timeix--time stampe for file; r_timxix--time stamp for the interested
     %intreval
     r_ix=start_ti:step_ti:stop_ti;
     pro_data(:,:)=NaN;
     for ix=1:length(timeix)
         %assign every element of data to the correct position
         pro_data(find(r_ix==timeix(ix)),:)=dataLevel(ix,:);
     end
   end
    
    %convert to duty cycle
    [cut_point]=Recursive_oneside_hypthesis_testing(dataLevel, 100);
    thres=cut_point(end);
    cs=pro_data>thres;
    scr=sum(cs,2)/length(cs(1,:));
    %set NaN element 
    scr(find(isnan(pro_data(:,1))))=NaN;
    
    %save scr data
    data_sp(:,n)=scr(:,1); %make sure the two vector matches
    locateion_sp(1,n)=deviceInfo.Longitude;
    locateion_sp(2,n)=deviceInfo.Latitude;
end

display('finish spatial data load!!!');

%% re-save loaded data
save_path='D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/%s';
save(sprintf(save_path, 'Spdata_1730_1740.mat'),'data_sp');
save(sprintf(save_path,'Sptime_1730_1740.mat'),'dateStamp');
save(sprintf(save_path,'Splocation_1730_1740.mat'),'locateion_sp');