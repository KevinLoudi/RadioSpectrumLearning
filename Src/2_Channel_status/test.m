% Created on Mon 26th Feb  13:25:45 2017
% Propose: Test script for part 2. channel status deciding
% Enviroment: Matlab 2015b
% @auththor: kevin

load D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1720_1760.mat;
%% thresholding
local_data=dataLevel(1:100,1:100);
ThresInfo=DoubleThresholding(local_data);

%% plot histogram of the two classified group and fit into normal distributions

clf;
%plot sample of two classes separtely
s1=local_data(~ThresInfo.Mask);
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
s2=local_data(ThresInfo.Mask);
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


%% double-threshold, separate into three region
thres_1=x1.max;
thres_2=x2.max


