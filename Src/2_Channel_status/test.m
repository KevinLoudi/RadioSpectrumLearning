% Created on Mon 26th Feb  13:25:45 2017
% Propose: Test script for part 2. channel status deciding
% Enviroment: Matlab 2015b
% @auththor: kevin

clear; clc; close all;
load D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1740_1760.mat;
display('successfully load data set..');
%% thresholding
startix=(1740-1740)/0.025+1;
stopix=(1750-1740)/0.025+1;
local_data=dataLevel(1:1000,startix:stopix);
%normalized data to [0,255] for sake of image show
nor_data=Normalize_matrix(local_data,0,255);
figure(1); imagesc(nor_data); title('original data');
%ThresInfo=DoubleThresholding(nor_data);
Cutpoint=Recursive_oneside_hypthesis_testing(nor_data,20);
%imagesc(ThresInfo.Mask); %show the thresholding results
%% classify
mask=nor_data>Cutpoint(end);
imagesc(mask);
save('ChannelStatusDataset/RandomCS_1740_1750.mat','mask')
imagesc(mask);  title('Channel status'); xlabel('Frequency'); ylabel('Time slot');

%% plot histogram of the two classified group and fit into normal distributions
%plot sample of two classes separtely
figure(2);
s1=nor_data(~ThresInfo.Mask);
h1=histogram(s1,'Normalization','probability');
hold on;

%fit distribution and plot pdf
pd1 = fitdist(s1,'Normal');
x1.min = 0;%pd1.mu-3.*(pd1.sigma);
x1.max = ThresInfo.Level;%pd1.mu+3*pd1.sigma;
%plot fitted distribution
x1.value=x1.min:((x1.max-x1.min)/100):x1.max
y1.value=pdf(pd1,x1.value);
hold on;
plot(x1.value, y1.value, 'LineWidth',2);

%plot sample of two classes separtely
s2=nor_data(ThresInfo.Mask);
h2=histogram(s2,'Normalization','probability');

%fit distribution and plot pdf
pd2 = fitdist(s2,'Normal');
x2.min = ThresInfo.Level;%pd2.mu-3.*(pd2.sigma);
x2.max = 255;%pd2.mu+3*pd2.sigma;
%plot fitted distribution
x2.value=x2.min:((x2.max-x2.min)/100):x2.max
y2.value=pdf(pd2,x2.value);
hold on;
plot(x2.value, y2.value, 'LineWidth',2);
title('classify analysis');

%% double-threshold, separate into three region
thres_1=pd1.mu+2*pd1.sigma;
thres_2=pd2.mu-2*pd2.sigma;
thri_results=zeros(size(nor_data));
[row,col]=size(thri_results);
for i=1:row
    for j=1:col
        if nor_data(i,j)<thres_1
            thri_results(i,j)=0;
        elseif nor_data(i,j)>thres_2
            thri_results(i,j)=2;
        else
            thri_results(i,j)=1; %uncertainty resion
        end
    end
end

region_noise=thri_results==0;
region_signal=thri_results==2;
region_uncertainty=thri_results==1;
figure(3); imagesc(region_noise); title('noise region');
figure(4); imagesc(region_signal); title('signal region');
figure(5); imagesc(region_uncertainty); title('uncertainty region');

%% test simulation data
clear; clc; close all;
load D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SimulationDataset\\Simset.mat;
display('sucessfully load sim data...');
%% plot data
figure(1);
plot(1:length(y),y);
hold on;
plot(1:length(rcs), (rcs>0)*100, 'LineWidth',2);

%% thresholding analysis
nor_y=Normalize_matrix(y,0,255);
thre_info=DoubleThresholding(y);
axis([0 120 0 0.2]);
% hold on;
% text(thre_info.LowThreshold,0,int2str(thre_info.LowThreshold));
% hold on;
% text(thre_info.HighThreshold,0,int2str(thre_info.HighThreshold));

%% classify
figure(2)
res=Classification(nor_y,thre_info);
err=sum(abs(res-rcs))/length(rcs); %  96.2%



