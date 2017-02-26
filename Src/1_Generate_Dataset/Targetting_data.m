% ****************************************************
% Note: Read *.argus file through a given filepath 
%           FreqInfo: startfreq, stopfreq, stepfreq
%           DeviceInfo: automaticly set 
%           TargetInfo: target-startfreq, target-stopfreq,target-stepfreq
%  Author: Kevin
%  Date: 24th Feb, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

function [Data,DeviceInfo]=Targetting_data(file_path,FreqInfo,TargetInfo)
 %device information and data
 Data.time='';
 Data.level=0;
 DeviceInfo.DeviceId=0;
 DeviceInfo.Longitude=0;
 DeviceInfo.Latitude=0;
 DeviceInfo.Status=0;
 DeviceInfo.FileNum=0;
 
 [DeviceInfo,Data]=ArgusReader_loc(file_path,FreqInfo,DeviceInfo,Data)

 %target part index
 low_t=(TargetInfo.StartF-FreqInfo.StartF)/FreqInfo.StepF;
 high_t=(TargetInfo.StopF-FreqInfo.StartF)/FreqInfo.StepF;
 
 Data.level=Data.level(low_t:high_t);
end

%Read all data sets from an argus file
function [Info,Data]=ArgusReader_loc(Path,FreqInfo,Info,Data)
 len = (FreqInfo.StopF-FreqInfo.StartF)/FreqInfo.StepF;
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
 jump_distance = 36+4;
 fseek(fid,jump_distance,'bof');
 Data.level(Info.FileNum, 1:len)= fread(fid,len,'integer*2=>int16',6);
 fclose(fid);
end