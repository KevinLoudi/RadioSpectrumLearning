%% 
% clear; clc; close all;
% len=1e5;
% [s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len);
% seq_str=strjoin(rcs,'');
% seq_mean=sum(rcs)/length(rcs);
% 
%% Preparation
alp = alphabet(seq_str);
ALGS = {'LZms', 'LZ78', 'PPMC', 'DCTW', 'PST'};;
params.ab_size = size(alp);
params.d = 10;  %maxmium order of VMM order %D
params.m = 2; %only for LZ-MS
params.s = 8; %only for LZ-MS
params.pMin = 0.006; %only for PST
params.alpha= 0; %only for PST
params.gamma = 0.0006; %only for PST
params.r = 1.05; %only for PST
params.vmmOrder = params.d; 

%split the dataset into five part
split_point=[1, floor(length(seq_str)/5),floor(2*length(seq_str)/5),floor(3*length(seq_str)/5),...
    floor(4*length(seq_str)/5), length(seq_str)];

datasets=cell(1,5);
for i=1:5
    datasets{i}=seq_str(split_point(i):split_point(i+1));
%     display(sprintf('Dataset %s',num2str(i)));
%     size(datasets{i});
end
train_sets=strcat(datasets{1},datasets{2},datasets{4},datasets{5});
test_sets=datasets{3};

pre_steps=20;
rand_sel=floor((length(test_sets)-pre_steps)*rand());
context=test_sets((rand_sel-pre_steps+1):rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+pre_steps));
pre_states=zeros(1,pre_steps);


% %% Model match level
% %algorithm 
% disp('---------------------------------------------------');
% disp('working with AB={0,1}');
% disp('---------------------------------------------------');
% disp(' ');
% 
% % 3. run each of the VMM algorithms
% for i=1:length(ALGS),
%     disp(sprintf('Working with %s', ALGS{i} ));
%     disp('--------')
%     %create a VMM through training 
%     t=cputime;
%     jVmm = vmm_create(map(alp, train_sets),  ALGS{i}, params);
%     
%     disp(sprintf('Pr(0 | 0001) = %f', vmm_getPr(jVmm, map(alp,'0'), map(alp,'0001'))));
%     e=double(cputime-t),
%     %disp(sprintf('Pr(1 | 1110) = %f', vmm_getPr(jVmm, map(alp,'1'), map(alp,'1110'))));
%     % calculates the length in bits of the  "compressed" representation of
%     % seq.  -log[ Pr ( seq | jVmm) ]
%     disp(sprintf('-lg(Pr(tar))=%f', vmm_logEval(jVmm,map(alp, context))));
%     disp('--------')
%     disp(' ');
% end

%% Prediction ability


%% 
display('Test prediction ability....');

for ii=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{ii} ));
    disp('--------')
  for i=1:pre_steps
    jVmm = vmm_create(map(alp, train_sets),  ALGS{ii}, params);
    p_0=vmm_getPr(jVmm, map(alp,'0'), map(alp,context));
    p_1=vmm_getPr(jVmm, map(alp,'1'), map(alp,context));
    if p_0>p_1
        pre_states(i)=0;
    else
        pre_states(i)=1;
    end
    %ipdate context
    context(1)=[];
    context=strcat(context, rcs(i));
  end
  display(sprintf('Prediction results of %s', ALGS{ii}));
  pre_cs=strjoin(pre_states,'')
  err=(pre_cs-rcs);
  display(sprintf('Error rate of %s: ', ALGS{ii}, num2str(sum(abs(err))/length(err))));
  
end


