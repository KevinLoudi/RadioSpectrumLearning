% Created on WED Mar 8th 13:25:45 2017
% Propose: Thresholding with a sliding window
% Enviroment: Matlab 2015b
% @auththor: kevin

function [cs,thres]=SlidingThresholding(data,freqix,freqwidth)

width=freqwidth/0.025; %40
startix=1; ix=startix;
stopix=startix+length(freqix);
cs=zeros(size(data));
thres=[];

while ix<=stopix-width
     local_data=data(:,ix:(ix+width));
     [cut_point]=Recursive_oneside_hypthesis_testing(local_data,20);
     local_thres=cut_point(end);
     thres=[thres,local_thres];
     cs(:,ix:(ix+width))=data(:,ix:(ix+width))>local_thres; %local thresholding
     ix=ix+width; %move to the next region
end

local_data=data(:,ix:end);
[cut_point]=Recursive_oneside_hypthesis_testing(local_data,20);
local_thres=cut_point(end);
thres=[thres,local_thres];
cs(:,ix:end)=data(:,ix:end)>local_thres; %local thresholding
end



