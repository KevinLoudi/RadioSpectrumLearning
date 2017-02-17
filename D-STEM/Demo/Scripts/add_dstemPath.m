actualPath  =  fileparts(which(mfilename));
if ispc
    idx = strfind(actualPath,'\');
    dstemPath = [actualPath(1:idx(end-1)), 'Src'];
else
    idx = strfind(actualPath,'/');
    dstemPath = [actualPath(1:idx(end)), 'Src'];
end
addpath(dstemPath);