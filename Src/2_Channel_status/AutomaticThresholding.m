% Created on Mon 27th Feb  13:25:45 2017
% Propose: Automatic thresholding data into three regions
% Enviroment: Matlab 2015b
% @auththor: kevin
% ThresInfo.LevelArr, Level, Mask

%% 
clear; clc; close all;
load D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_80_110.mat;
display('successfully load data set..');
%dataLevel=Normalize_matrix(dataLevel,0,255);

%%
StartF=80; %MHz
StopF=110;

win_f=40; %freqency width 1MHz, compare to 0.25MHz bandwidth, 40
win_t=60; %time length 4hours, take 2400 samples into count

[f_len,t_len]=size(dataLevel);

for f=(win_f+1):f_len
    for t=(win_t+1):t_len
        win=dataLevel((f-win_f):(f-1), (t-win_t):(t-1));
        display('inital a observe window!!!');
        %calculate threshold
        [thres_info]=DoubleThresholding(win);
    end
end




