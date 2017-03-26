% Created on Mon 25th Feb 13:25:45 2017
% Propose: Load mutiply sensors' Spectrum Data from Binary Files
% Enviroment: Matlab 2015b
% @auththor: kevin

function LoadMultiSensorsSpectrum()
%%Enabled when file read is needed
 %File Folder
 recalculate_all=0; %skip the loaded data set 
 display('Load multiply sensors'' data of a given frequency range.');
 display('Step1: set file path and data load parameters....');
 SensorIds={'02','03','04','06','07','08','09','10','11','12','13','15','17','18',...
     '20','21','22','23','25','26','27','28','30','31','33','34','35','36','37','38','39',...
     '40','42','45','46','47','48','49','51','52','53','54','55','56','58','59','60','61',...
     '62','63','64','65','66','68','69','71','72','73','76','77','78','79','80','81','82',...
     '83','84','85','87','88','89','90','92'};

 %with trobule: '19' '50'
 
 Path ='D:\\Code\\Data\\60-137\\%s\\';
 Def_path=Path;
 
display('Step2: read file and load data....');
 for i=1:length(SensorIds)
   display('=====================');
   display('deal with a sensor....');
 
   Path=Def_path;
   Path=strcat(Path, SensorIds{i});
   
   %Path ='D:\\Code\\Backup\\Matlab\\SpectrumLearning\\Data\\OriginData\\1710-1960\\%s\\02';
    %Critical Parameters
    StartF=88000;%kHz
    StopF=108000;%
    StepF = 25;
    
   check_name=sprintf('MultiSensorsDataset\\Dataset_%s_%s_%s.mat',int2str(SensorIds{i}), ...
      int2str(StartF),int2str(StopF));
   if exist(check_name) && ~recalculate_all
       continue; %skip this loop and go to see the next device
   end
    
    DayArray = {'20151216','20151217','20151218','20151219','20151220','20151221','20151222'};
    
    %Read data day by day if the "Path" exist locally
    
    [MultiLevel] = MultiDaySpectrumReader(Path,DayArray,StartF,StopF);
 
    display('shape data and cut off....');
    MultiLevel.Info.StartFreq = StartF;
    MultiLevel.Info.StopFreq = StopF;
    MultiLevel.Info.StepFreq = StepF;
    cur_time = fix(clock);
    time_str = sprintf('%.4d-%.2d-%.2d:%.2d:%.2d:%.2d:%.2d',cur_time(1),cur_time(2),cur_time(3),cur_time(4),cur_time(5),cur_time(6));
    MultiLevel.Info.BuildTime = time_str;
    filename = sprintf('MultiLevel_%s_%s.mat', num2str(StartF),num2str(StopF));
    save(filename,'MultiLevel');

    %Write MultiLevel Data and Time Stamps into CSV file
    %  filename='MultiLevel_60_137.mat';
    display('save data in a unit pack....');
    CSVFormatWriter(filename,7);
 end
 disp('Step3: successfully output dataset!!!');
end

function CSVFormatWriter(FileName,days)
  load(FileName);
  dateStamp=[];
  dataLevel=[];
  for d=7
      timeStamparr=MultiLevel.ByDay{d}.time;
      levelDataarr=MultiLevel.ByDay{d}.level;
      %add level data from the new day in matrix
      dataLevel=[dataLevel;levelDataarr];
      %deal with missing data from the vary beginning
      for i = 1:length(timeStamparr)
        tmpTime=timeStamparr(i,1:12);
        tmpDataVector=[str2num(tmpTime(1:4)),str2num(tmpTime(5:6)),str2num(tmpTime(7:8)),...
            str2num(tmpTime(9:10)),str2num(tmpTime(11:12)),0];
        %add time stamp of the new day/timesolts
        dateStamp=[dateStamp; datestr(tmpDataVector,'yyyy-mm-dd HH:MM:SS')];
        %dateStamp=[dateStamp; datestr(tmpDataVector,'yyyy-mm-dd HH:MM:SS')];
     end
  end
  %save('timestamp.mat','dateStamp');
  %save('level.mat','dataLevel');
  %name protocol: Dataset_deviceid_startfreq_stopfreq
  dataSetname=sprintf('MultiSensorsDataset\\Dataset_%s_%s_%s.mat',int2str(MultiLevel.Info.DeviceId), ...
      int2str(MultiLevel.Info.StartFreq),...
      int2str(MultiLevel.Info.StopFreq));
  deviceInfo=MultiLevel.Info;
  
%   timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
  
%   start_ti=datetime('2015-12-15 00:00:00','InputFormat','yyyy-MM-dd HH:mm:SS');
%   step_ti=minutes(5); %5 minutes duration between 2 time slots
%   stop_ti=datetime('2015-12-22 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS');
%   freqix=MultiLevel.Info.StartFreq:MultiLevel.Info.StopFreq:MultiLevel.Info.StopFreq;
  
%   exp_ti_slots=(stop_ti-start_ti)/step_ti+1; %expected time slots
%   pro_data=0;
%   dataLevel=sum(dataLevel,2)/length(dataLevel(1,:));
%   %deal with missing data
%   if exp_ti_slots==length(dataLevel(:,1))
%      display('no missed data');
%      pro_data=dataLevel;
%  %unequal to the therotical value, must miss some data
%  else
%      display('need to process missed data!!!');
%      pro_data=zeros(exp_ti_slots,length(freqix)); %prepare to process the missing data
%      ti=start_ti; ixt=1; %index for time slots
%      ix=1; %index for the row of dataset
%      
%      %timeix--time stampe for file; r_timxix--time stamp for the interested
%      %intreval
%      r_ix=start_ti:step_ti:stop_ti;
%      pro_data(:,:)=NaN;
%      for ix=1:length(timeix)
%          %assign every element of data to the correct position
%          pro_data(find(r_ix==timeix(ix)),:)=dataLevel(ix,:);
%      end
%   end
% 
%   %renew data and timestamp
%   dataLevel=pro_data;
%   dateStamp=datestr(r_ix,'yyyy-mm-dd HH:MM:SS'); 
%   clear pro_data;
 
  save(dataSetname,'deviceInfo','dateStamp','dataLevel');
  display('finish work on a sensor..');
  display('=====================');
end

%Read Out Spectrum Data from several days (one station)
%Author: Zhu Gengyu
%Date: 2016/7/18
%Path format: 'Y:\\20-3000\\%s\\03' DayArray should be string
%Level.Info: Device, Lon, Lat, Status, FileNum  Level.Data: Time[time,stringlen], level[time,freq]

function [MultiLevel] = MultiDaySpectrumReader(Path,DayArray,StartFreq,StopFreq)
Maxsize = 1000; 
%index for the time slot num
datasize = 1;
%declear Level container
%Level.Info=0;  %Device, Lon, Lat 
Data.time='';
Data.level=0;
Info.DeviceId=0;
Info.Longitude=0;
Info.Latitude=0;
Info.Status=0;
Info.FileNum=0;
for i = 1:length(DayArray) 
  Day = DayArray{i};  
  %get the file path for a specific day
  DayPath = sprintf(Path,Day);
  %find all files in the path
  dirinfo = dir(DayPath);
  index = 1;
  %Loop for all files within the path
  for k = 1:length(dirinfo)
     %point to this-file
     thisdir = dirinfo(k).name;
     %get the path for this-file
     filename = [DayPath,'\',thisdir];
     %check if this-file really exist
     if exist(filename,'file') == 2
        %Read out data from a *.argus file
        [Info,Data]=ArgusReader(filename,StartFreq,StopFreq,Info,Data);
        %tmpLevel.SampleTime=floor(str2num(thisdir(1:12)));
        Data.time(Info.FileNum,1:12)=thisdir(1:12);
        %Count the accessed data
        datasize = datasize + 1;
        index = index + 1;
        if(index > Maxsize)
            %abortion if the data set is too large
            delete Info, Data ;
            break;
        end
     end
  end 
  %Combine the newly read data and existed data in Level
  MultiLevel.ByDay{i}=Data;
end
%presever other information
MultiLevel.Info=Info;
end

%Read all data sets from an argus file
function [Info,Data]=ArgusReader(Path,StartF,StopF,Info,Data)
 
 len = (StopF-StartF)/25+1;
 fid = fopen(Path);
 jump_distance = 0;
 fseek(fid,jump_distance,'bof');
 Info.FileNum=Info.FileNum+1;
 if  Info.Status==0
   Info.DeviceId=fread(fid,1,'integer*4=>int32');
   Info.Longitude= fread(fid,1,'float=>float');
   Info.Latitude=fread(fid,1,'float=>float');
   Info.Status=1; %aleard assign device information
 end

%  jump_distance = 36+0; %+0: min level; +2: mean level; +4: max level
%  fseek(fid,jump_distance,'bof');
%  LevelData.MinLevel = fread(fid,len,'integer*2=>int16',6);
%  jump_distance = 36+2;
%  fseek(fid,jump_distance,'bof');
%  LevelData.MeanLevel= fread(fid,len,'integer*2=>int16',6);
 jump_distance = 36+4;
 fseek(fid,jump_distance,'bof');
 %LevelData.MaxLevel= fread(fid,len,'integer*2=>int16',6);
 Data.level(Info.FileNum, 1:len)= fread(fid,len,'integer*2=>int16',6);
 fclose(fid);
end