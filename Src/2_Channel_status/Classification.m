% Created on Mon 26th Feb  13:25:45 2017
% Propose: do classfication according to the double threshold
% Enviroment: Matlab 2015b
% @auththor: kevin
% ThresInfo.LevelArr, Level, Mask, LowThreshold, HighThreshold

function [Res]=Classification(Data, Threshold)
    if nargin<2
        error('Err: Lack input parameters!!!');
    elseif max(Data(:))<Threshold.LowThreshold || min(Data(:))>Threshold.HighThreshold
        error('Err: Threshold do not match the data set!!!');
    elseif Threshold.LowThreshold>Threshold.HighThreshold
        warning('Warning: threshold mismatch!!!');
        t2=Threshold.LowThreshold;
        t1=Threshold.HighThreshold;
    else
        t1=Threshold.LowThreshold;
        t2=Threshold.HighThreshold;
    end

    [row,col]=size(Data);
    Res=zeros(size(Data));
    %do classification
    for i=1:row
     for j=1:col
        if Data(i,j)<t1
            Res(i,j)=0;
        elseif Data(i,j)>t2
            Res(i,j)=2;
        else
            Res(i,j)=1; %uncertainty resion
        end
       end
    end

end