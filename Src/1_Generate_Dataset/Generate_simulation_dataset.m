% Created on Mon 27th Feb  13:25:45 2017
% Propose: Generate simulation data to check the performance of 
%   detection and thresholding ability
% Enviroment: Matlab 2015b
% @auththor: kevin

%clear; thresh=Generate_signal();
function [y,rcs]=Generate_simulation_dataset()
    clear;
    %generate signal
    len=10000;
    t=1:len;
    
    %noise with mean and std
    n=Uniform_noise(5, 2, len);
    %PU signal with mean, std and expected visit times
    [s1,rcs1]=Singal_source(35,10,len,8);
    [s2,rcs2]=Singal_source(36,10,len,5);
    %[s3,rcs3]=Singal_source(24,10,len,2);
    
    %SU received signal
    y=s1+s2+n;
    rcs=rcs1+rcs2;
    
    %plot(1:len, y);
end

function n=Uniform_noise(mean, std, len)
    n=random('Normal', mean, std, 1, len);
end

%simulate a singal source
%s: generated signal, x: source presence state
%mean, std, len: mean, standard variance and length of signal
%visit_time: occur times of the signal
function [s,x]=Singal_source(mean, std, len, visit_time)
    t=1:len;
    
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

    %generate signal when it present
    s=random('Normal', mean, std, 1, len).*x;
end
