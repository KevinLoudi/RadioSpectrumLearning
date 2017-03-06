% Created on Mon 26th Feb  13:25:45 2017
% Propose: Test script for part 3. time series analysis
% Enviroment: Matlab 2015b
% @auththor: kevin

%% Load data
clear; clc; close all;
load D:\\Code\\WorkSpace\\ThesisCode\\Src\\2_Channel_status\\ChannelStatusDataset\\RandomCS_1740_1750.mat;
display('successfully load channel status data set..');
figure(1);
imagesc(mask);  title('Channel status'); xlabel('Frequency'); ylabel('Time slot');

%% Shape data--define sequence and its alphabet
seq=int16(mask(:,1)');
%assemble a string
seq_str='';
for i=1:length(seq)
    seq_str=strcat(seq_str,int2str(seq(i)));
end
tar='0';
len=length(seq);
% figure(2);
% plot(1:len,seq);  axis([1 len+1 0 2]);
ab = alphabet(seq_str);
split_point=980;
seq_train=seq_str(1:split_point);
seq_test=seq_str((split_point+1):end);

%% Test several algorithms
%ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'BinaryCTW', 'PST'};
ALGS = {'PPMC', 'BinaryCTW'};
params.ab_size = size(ab);
params.d = 10;  %two state
params.m = 2;
params.s = 8;
params.pMin = 0.006;
params.alpha= 0;
params.gamma = 0.0006;
params.r = 1.05;
params.vmmOrder = params.d;

disp('---------------------------------------------------');
disp('working with AB={0,1}');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    %create a VMM through training 
    jVmm = vmm_create(map(ab, seq_train),  ALGS{i}, params);
    disp(sprintf('Pr(0 | 00) = %f', vmm_getPr(jVmm, map(ab,'0'), map(ab,'00'))));
    disp(sprintf('Pr(1 | 00) = %f', vmm_getPr(jVmm, map(ab,'1'), map(ab,'00'))));
    %averaged log-loss
    disp(sprintf('-lg(Pr(tar))=%f', vmm_logEval(jVmm,map(ab, seq_test))));
    disp('--------')
    disp(' ');
end


%% Load data for DC calculation
cut_point=(1745-1740)/0.025;
mask=mask(:,1:cut_point);
figure(3)
imagesc(mask);  title('Channel status'); xlabel('Frequency'); ylabel('Time slot');

%% Calculate DC data
[tn,fn]=size(mask); dc=zeros(tn,1);
for i=1:tn
    dc(i,1)=sum(mask(i,:))/fn;
end
figure(4)
plot(1:tn,dc); title('Duty cycle');

csvwrite('TimeSeriesDataset/dc_2_1740_1750.csv',dc);


