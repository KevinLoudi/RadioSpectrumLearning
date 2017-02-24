% ****************************************************
% Note: Generate a dataset for further usage
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

file_path='D:\\Code\\Backup\\Matlab\\SpectrumLearning\\Data\\OriginData\\1710-1960\\20151216\\02\\201512151900.argus';
 
FreqInfo.StartF= 1710;%60;
FreqInfo.StopF=1960;%137;
FreqInfo.StepF = 0.025;

TargetInfo.StartF= 1720;%60;
TargetInfo.StopF=1740;%137;
TargetInfo.StepF = 0.025;

[Data,DeviceInfo]=Targetting_data(file_path,FreqInfo,TargetInfo);


%Critical Parameters
 StartF= 1710;%60;
 StopF=1760;%137;
 StepF = 0.025;

 %device information
 Data.time='';
 Data.level=0;
 Info.DeviceId=0;
 Info.Longitude=0;
 Info.Latitude=0;
 Info.Status=0;
 Info.FileNum=0;
 
 [Info,Data]=ArgusReader(file_path,StartF,StopF,Info,Data);