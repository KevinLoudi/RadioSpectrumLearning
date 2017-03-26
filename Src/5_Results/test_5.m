%% Prediction ability influenced by order
clear; clc; close all;
len=1e5;
[s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len,10,0.5,4);
seq_str=strjoin(rcs,'');
seq_mean=sum(rcs)/length(rcs);

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
end
train_sets=strcat(datasets{1},datasets{2},datasets{4},datasets{5});
test_sets=datasets{3};

pre_steps=100; %prediction length
train_steps=20000; %80000

rand_sel=floor((length(test_sets)-pre_steps)*rand());
context=test_sets(rand_sel-pre_steps:rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+pre_steps));
pre_states=zeros(1,pre_steps);

D_arr=[5 10 20 30 40];
train_steps_arr=100:200:80000;

for j=1:length(D_arr)
    params.d =D_arr(j)
    
for f=1:length(train_steps_arr)
    train_steps=train_steps_arr(f);
    t_train_sets=train_sets(1:train_steps);
   
for ii=3%1:length(ALGS),
    tstart_1=tic;
    jVmm = vmm_create(map(alp, t_train_sets),  ALGS{ii}, params);
    tend_1=toc(tstart_1);
%     disp(sprintf('Trainning time %s', num2str(tend_1) ));
    train_time_arr(f,j)=tend_1;
    
    tstart_2=tic;
    for i=1:pre_steps
    
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
     tend_2=toc(tstart_2);
%      disp(sprintf('Prediction time %s', num2str(tend_2) ));
     predict_time_arr(f,j)=tend_2;
%   display(sprintf('Prediction results of %s', ALGS{ii}));
  pre_cs=strjoin(pre_states,'');
  err=(pre_cs-rcs);
  err_rate=sum(abs(err))/length(err);
  err_arr(f,j)=err_rate;
%   display(sprintf('Error rate of %s: ', ALGS{ii}, num2str(err_rate)));
  
end
end
end
%% 
D_arr=[5 10 20 30 40];
figure(1);
subplot(2,1,1);plot(train_steps_arr,err_arr,'LineWidth',1.0);
legend('D = 5','D = 10','D = 20','D = 30','D = 40', 'FontSize',12);
ylabel('Ô¤²âÎó²î', 'FontSize',12);
subplot(2,1,2);plot(train_steps_arr,train_time_arr,'LineWidth',1.0);
legend('D = 5','D = 10','D = 20','D = 30','D = 40', 'FontSize',12);
xlabel('ÑµÁ·²½Êý', 'FontSize',12); ylabel('ÑµÁ·Ê±³¤/sec', 'FontSize',12);
path='D:/doc/PapaerLibrary/Figures/vmm_order';
% print(path,'-dpng','-r500');