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
seq=mask(:,1);
tar='10001';
len=length(seq);
figure(2);
plot(1:len,seq);  axis([1 len+1 0 2]);
ab = alphabet('10001');

%% Test several algorithms
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'BinaryCTW', 'PST'};
params.ab_size = size(ab);
params.d = 2;  %two state
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
    jVmm = vmm_create(map(ab, seq),  ALGS{i}, params);

    disp(sprintf('Pr(0 | 000) = %f', vmm_getPr(jVmm, map(ab,'0'), map(ab,'000'))));
    disp(sprintf('Pr(1 | 000) = %f', vmm_getPr(jVmm, map(ab,'1'), map(ab,'000'))));;

    disp(sprintf('-lg(Pr(tar))=%f', vmm_logEval(jVmm,map(ab, tar))));
    disp('--------')
    disp(' ');
end

%% Demo of the toolbox

% A simple hands-on tutorial 
% Browse the code to get the basics of how-to utilize this VMM tool

% 1. defining sequence and its alphabet
seq = 'abracadabra';
ab = alphabet('abracadabra');

% 2. testing all algs; param values ~match "best" values for text data (see
% tbl.8 in VMM paper)
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'BinaryCTW', 'PST'};
params.ab_size = size(ab);
params.d = 5;
params.m = 2;
params.s = 8;
params.pMin = 0.006;
params.alpha= 0;
params.gamma = 0.0006;
params.r = 1.05;
params.vmmOrder = params.d;

% use AB with size = 5
disp('---------------------------------------------------');
disp('working with AB={a, b, c, d, r }');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    jVmm = vmm_create(map(ab, seq),  ALGS{i}, params);

    disp(sprintf('Pr(a | br) = %f', vmm_getPr(jVmm, map(ab,'a'), map(ab,'br'))));
    disp(sprintf('Pr(b | br) = %f', vmm_getPr(jVmm, map(ab,'b'), map(ab,'br'))));
    disp(sprintf('Pr(c | br) = %f', vmm_getPr(jVmm, map(ab,'c'), map(ab,'br'))));
    disp(sprintf('Pr(d | br) = %f', vmm_getPr(jVmm, map(ab,'d'), map(ab,'br'))));
    disp(sprintf('Pr(r | br) = %f', vmm_getPr(jVmm, map(ab,'r'), map(ab,'br'))));

    disp(sprintf('-lg(Pr(abracadabra))=%f', vmm_logEval(jVmm,map(ab, 'abracadabra'))));
    disp('--------')
    disp(' ');
end

% 4. repeat the same scenario - this time with ascii AB (size(AB)=127)

disp(' ');
disp(' ');
disp('---------------------------------------------------');
disp('working with ascii AB (|AB|=127)');
disp('---------------------------------------------------');
disp(' ');
params.ab_size = 127;


for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    jVmm = vmm_create( seq,  ALGS{i}, params);

    disp(sprintf('Pr(a | br) = %f', vmm_getPr(jVmm, 'a', 'br')));
    disp(sprintf('Pr(b | br) = %f', vmm_getPr(jVmm, 'b', 'br')));
    disp(sprintf('Pr(c | br) = %f', vmm_getPr(jVmm, 'c', 'br')));
    disp(sprintf('Pr(d | br) = %f', vmm_getPr(jVmm, 'd', 'br')));
    disp(sprintf('Pr(r | br) = %f', vmm_getPr(jVmm, 'r', 'br')));

    disp(sprintf('-lg(Pr(abracadabra))=%f', vmm_logEval(jVmm, 'abracadabra')));
    disp('--------')
    disp(' ');
end
