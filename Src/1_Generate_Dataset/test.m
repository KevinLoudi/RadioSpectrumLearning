% Created on Mon 27th Feb  13:25:45 2017
% Propose: Test code for generating dataset
% Enviroment: Matlab 2015b
% @auththor: kevin

%% generate simulation data
clear;clc;close all;
[y,rcs]=Generate_simulation_dataset();
%% plot data
figure(1);
plot(1:length(y),y);
hold on;
plot(1:length(rcs), (rcs>0)*100, 'LineWidth',2);
save('SimulationDataset/Simset.mat','y','rcs');
