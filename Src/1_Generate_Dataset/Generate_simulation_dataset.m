% Created on Mon 27th Feb  13:25:45 2017
% Propose: Generate simulation data to check the performance of 
%   detection and thresholding ability
% Enviroment: Matlab 2015b
% @auththor: kevin

%clear; thresh=Generate_signal();  SourceInfo is a cell
function [y,cs]=Generate_simulation_dataset(SourceInfo,NoiseInfo,len)
    %generate signal
    t=1:len;
    
    %noise with mean and std
    n=Uniform_noise(NoiseInfo.mean, NoiseInfo.std, len);
    y=0; cs=0;
    %PU signal with mean, std and expected visit times
    for i=1:length(SourceInfo)
        [s,rcs]=Singal_source(SourceInfo{i}.mean,SourceInfo{i}.std,len,SourceInfo{i}.visit_time);
        y=y+s;
        cs=cs+rcs;
    end
    %plus envirmental noise %SU received signal
    y=y+n;
    
    plot(1:len, y);
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
        a=1; b=max_possible_gap(i)/2; %random number boundry
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
