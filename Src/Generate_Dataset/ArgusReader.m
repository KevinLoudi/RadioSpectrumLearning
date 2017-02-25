
%Read all data sets from an argus file
function [Info,Data]=ArgusReader(Path,StartF,StopF,Info,Data)
 len = (StopF-StartF)/0.025+1;
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