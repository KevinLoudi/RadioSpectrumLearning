% ****************************************************
%  Note: Test script for the project
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

clear;clf;clc;
%generate signal
[data,rcs]=Generate_signal();
%select threshold
thresh=Threshold_selection_load(data);
cs=data>thresh(1);
%error rate 
rcs=rcs>0;
err_s=abs(rcs-cs);
err=sum(abs(rcs-cs))/length(data);

plot(1:length(data),err_s*100,'b');
hold on;
plot(1:length(data),data,'g');
% hold on;
% plot(1:length(data),(rcs>0)*100,'k');




