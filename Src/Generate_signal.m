% ****************************************************
% Note: This is a part of the final code aim at generating random signal
%    with random presence
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

clear;

%time
len=1000;
t=1:len;

%randomly set visit start time
visit_time=10;
pd=makedist('Normal','mu',0,'sigma',1);
%generate random number in range [a,b]
a=1; b=len;
start_time=floor((b-a)*rand(visit_time,1)+a);
start_time=start_time';
start_time=sort(start_time,'ascend');

%randomly set leave time 
max_possible_gap=zeros(1,visit_time);
%latest time for a signal to stay
latest_level_time=[start_time(2:end) , len];
for i=1:visit_time
    %calculate the maxiest possible signal stay time
    max_possible_gap(i)=latest_level_time(i)-start_time(i);
end
%randomly set stay time
stay_time=zeros(1,visit_time);
for i=1:visit_time
    a=1; b=max_possible_gap(i); %random number boundry
    %generate random number in range of [a,b]
    stay_time(i)=floor((b-a)*rand(1)+a);
end

%set signal presence state
x=zeros(1,len);
for i=1:visit_time
    %set the signal-presence span as '1'
    x(start_time(i):(start_time(i)+stay_time(i)))=x(start_time(i):(start_time(i)+stay_time(i)))+1;
end
%plot signal-presence status
plot(t,x);


