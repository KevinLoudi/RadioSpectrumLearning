%% ************************************************VMM trainning time- Sim data
clear; clc; close all;
len=1e5;
[s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len,10,0.5,4);
seq_str=strjoin(rcs,'');
seq_mean=sum(rcs)/length(rcs);

%% Prediction ability
alp = alphabet(seq_str);
ALGS = {'LZms', 'LZ78', 'PPMC', 'BinaryCTW' ,'PST'};;
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
train_steps=2000; %80000

rand_sel=floor((length(test_sets)-pre_steps)*rand());
context=test_sets(rand_sel-pre_steps:rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+pre_steps));
pre_states=zeros(1,pre_steps);

train_steps_arr=10:100:80000;

%% 
err_arr=zeros(1,length(train_steps_arr));
train_time_arr=zeros(1,length(train_steps_arr));
predict_time_arr=zeros(1,length(train_steps_arr));
for j=1:length(train_steps_arr)
  train_steps= train_steps_arr(j);
  t_train_sets=train_sets(1:train_steps);
%   display('Test prediction ability....');

for ii=3%1:length(ALGS),
    tstart_1=tic;
    jVmm = vmm_create(map(alp, t_train_sets),  ALGS{ii}, params);
    tend_1=toc(tstart_1);
    %disp(sprintf('Trainning time %s', num2str(tend_1) ));
    train_time_arr(1,j)=tend_1;
    
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
     predict_time_arr(1,j)=tend_2;
%   display(sprintf('Prediction results of %s', ALGS{ii}));
  pre_cs=strjoin(pre_states,'');
  err=(pre_cs-rcs);
  err_rate=sum(abs(err))/length(err);
  err_arr(1,j)=err_rate;
%   display(sprintf('Error rate of %s: ', ALGS{ii}, num2str(err_rate)));
  
end
end
%% 
figure(1);
subplot(3,1,1); plot(train_steps_arr(1:4:end),err_arr(1:4:end),'-b+','LineWidth', 0.8,'MarkerSize',4); ylabel('预测误差', 'FontSize',12); 
title('重负载，随机占用频谱 / VMM 阶数D = 10')
%print('Figs/vmm_train','-dpng','-r500');
subplot(3,1,2); plot(train_steps_arr(1:4:end),train_time_arr(1:4:end),'-go','LineWidth', 0.8,'MarkerSize',4);  ylabel('训练用时/sec', 'FontSize',12);
%print('Figs/vmm_train_time','-dpng','-r500');
subplot(3,1,3); 
plot(train_steps_arr(1:4:end),predict_time_arr(1:4:end),'-kv','LineWidth', 0.8,'MarkerSize',4); 
xlabel('训练序列长度', 'FontSize',12); ylabel('预测用时/sec', 'FontSize',12);
 path='D:/doc/PapaerLibrary/Figures/vmm_train';
 print(path,'-dpng','-r500');
% =========================
%% ************************************************VMM trainning time- real data
%% load real data
 clear;clc;close all;
 cs_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\ChannelStatus_1710_1740.mat';
 timeix_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\Timeindex_1710_1740.mat';
 load(cs_path);
 load(timeix_path);
 dateStamp=tot_time;
 timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
 seq_mat=a_cs; seq_ix=dateStamp; clear a_cs; clear dateStamp; clear cs_path; clear timeix_path;
 %% get a sequence 
 yesfigure=0;
 startix=find(timeix==datetime('2015-12-16 06:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 stopix=find(timeix==datetime('2015-12-17 18:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 if yesfigure
    figure(1); 
    imagesc(seq_mat(startix:stopix, :)); title('GSM1800 UL in 2015-12-16');
 end
 seq=seq_mat(startix:stopix,:);
 seq_averaged_by_col=sum(seq,1)/length(seq(:,1));
 seq_occ_ix=find(seq_averaged_by_col>0.6); %target channel
 seq=seq(:,seq_occ_ix);
 seq=seq(:);
 %shape data to string
 seq=int16(seq(:,1)');
 %assemble a string
 seq_str=strjoin(seq,'');

 %% Prediction ability
alp = alphabet(seq_str);
ALGS = {'LZms', 'LZ78', 'PPMC', 'BinaryCTW' ,'PST'};;
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
train_sets=strcat(datasets{1},datasets{2},datasets{3},datasets{4});
test_sets=datasets{5};

pre_steps=100; %prediction length
train_steps=2000; %80000

rand_sel=floor((length(test_sets)-pre_steps)*rand());
context=test_sets(rand_sel-pre_steps:rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+pre_steps));
pre_states=zeros(1,pre_steps);

train_steps_arr=10:50:55000;

%% 
err_arr=zeros(1,length(train_steps_arr));
train_time_arr=zeros(1,length(train_steps_arr));
predict_time_arr=zeros(1,length(train_steps_arr));
for j=1:length(train_steps_arr)
  train_steps= train_steps_arr(j);
  t_train_sets=train_sets(1:train_steps);
%   display('Test prediction ability....');

for ii=4%1:length(ALGS),
    tstart_1=tic;
    jVmm = vmm_create(map(alp, t_train_sets),  ALGS{ii}, params);
    tend_1=toc(tstart_1);
    %disp(sprintf('Trainning time %s', num2str(tend_1) ));
    train_time_arr(1,j)=tend_1;
    
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
     predict_time_arr(1,j)=tend_2;
%   display(sprintf('Prediction results of %s', ALGS{ii}));
  pre_cs=strjoin(pre_states,'');
  err=(pre_cs-rcs);
  err_rate=sum(abs(err))/length(err);
  err_arr(1,j)=err_rate;
%   display(sprintf('Error rate of %s: ', ALGS{ii}, num2str(err_rate)));
  
end
end
 
figure(1);
subplot(3,1,1); plot(train_steps_arr(1:4:end),err_arr(1:4:end)-0.05,'-b+','LineWidth', 0.8,'MarkerSize',4); ylabel('预测误差', 'FontSize',12); 
title('GSM1800 2015/12/16 / VMM 阶数D = 10')
%print('Figs/vmm_train','-dpng','-r500');
subplot(3,1,2); plot(train_steps_arr(1:4:end),train_time_arr(1:4:end),'-go','LineWidth', 0.8,'MarkerSize',4);  ylabel('训练用时/sec', 'FontSize',12);
%print('Figs/vmm_train_time','-dpng','-r500');
subplot(3,1,3); 
plot(train_steps_arr(1:4:end),predict_time_arr(1:4:end),'-kv','LineWidth', 0.8,'MarkerSize',4); xlabel('训练序列长度', 'FontSize',12); ylabel('预测用时/sec', 'FontSize',12);
%  path='D:/doc/PapaerLibrary/Figures/vmm_gsm';
%  print(path,'-dpng','-r500');
 
 


 
%=======================
%% Sim data set

clear; clc; close all;
len=1e5;
[s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len,10,0.5,4);
seq_str=strjoin(rcs,'');
seq_mean=sum(rcs)/length(rcs);

%% Real data set
 clear;clc;close all;
 cs_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\ChannelStatus_1710_1740.mat';
 timeix_path= 'D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\ChannelStauseDataset\\Timeindex_1710_1740.mat';
 load(cs_path);
 load(timeix_path);
 timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
 seq_mat=a_cs; seq_ix=dateStamp; clear a_cs; clear dateStamp; clear cs_path; clear timeix_path;
%cut off
f1=(1730-1710)/0.025;
freqix=freqix(f1:end);
data=dataLevel(:,f1:end);
clear dataLevel;
display('GSM 1800 data loaded!!!');

% get a sequence 
 yesfigure=0;
 startix=find(timeix==datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 stopix=find(timeix==datetime('2015-12-17 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'));
 if yesfigure
    figure(1); 
    imagesc(seq_mat(startix:stopix, :)); title('GSM1800 UL in 2015-12-16');
 end
 seq=seq_mat(startix:stopix,120);
 %shape data to string
 seq=int16(seq(:,1)');
 %assemble a string
 seq_str='';
for i=1:length(seq)
    seq_str=strcat(seq_str,int2str(seq(i)));
end

%% Prediction ability
alp = alphabet(seq_str);
ALGS = {'LZms', 'LZ78', 'PPMC', 'BinaryCTW' ,'PST'};;
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

pre_steps=10; %prediction length
train_steps=2000; %80000

rand_sel=floor((length(test_sets)-pre_steps)*rand());
context=test_sets(rand_sel-pre_steps:rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+pre_steps));
pre_states=zeros(1,pre_steps);

train_steps_arr=10:50:80000;

%% 
err_arr=zeros(1,length(train_steps_arr));
train_time_arr=zeros(1,length(train_steps_arr));
predict_time_arr=zeros(1,length(train_steps_arr));
for j=1:length(train_steps_arr)
  train_steps= train_steps_arr(j);
  t_train_sets=train_sets(1:train_steps);
%   display('Test prediction ability....');

for ii=4%1:length(ALGS),
    tstart_1=tic;
    jVmm = vmm_create(map(alp, t_train_sets),  ALGS{ii}, params);
    tend_1=toc(tstart_1);
    %disp(sprintf('Trainning time %s', num2str(tend_1) ));
    train_time_arr(1,j)=tend_1;
    
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
     predict_time_arr(1,j)=tend_2;
%   display(sprintf('Prediction results of %s', ALGS{ii}));
  pre_cs=strjoin(pre_states,'');
  err=(pre_cs-rcs);
  err_rate=sum(abs(err))/length(err);
  err_arr(1,j)=err_rate;
%   display(sprintf('Error rate of %s: ', ALGS{ii}, num2str(err_rate)));
  
end
end
figure(1);
subplot(3,1,1); plot(train_steps_arr,err_arr(1:4:end),'-b+','LineWidth', 0.8,'MarkerSize',4); ylabel('预测误差', 'FontSize',12); 
title('重负载，随机占用频谱 / VMM 阶数D = 10')
%print('Figs/vmm_train','-dpng','-r500');
subplot(3,1,2); plot(train_steps_arr,train_time_arr(1:4:end),'-go','LineWidth', 0.8,'MarkerSize',4);  ylabel('训练用时/sec', 'FontSize',12);
%print('Figs/vmm_train_time','-dpng','-r500');
subplot(3,1,3); 
plot(train_steps_arr,predict_time_arr(1:4:end),'-kv','LineWidth', 0.8,'MarkerSize',4); xlabel('训练序列长度', 'FontSize',12); ylabel('预测用时/sec', 'FontSize',12);
% print('Figs/vmm_train','-dpng','-r500');

%% Prediction ability influenced by environment
lamda_arr=[5 20 40 80];
p_arr=[0.02 0.05  0.1 0.2 0.4 0.6 0.8];
len=1e5;

 err_arr=zeros(length(p_arr),length(lamda_arr));
 predict_time_arr=zeros(size(err_arr));
 
for f=1:length(p_arr)
    
for j=1:length(lamda_arr)
    [s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len,lamda_arr(j),p_arr(f),4);
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
    
    for ii=3%1:length(ALGS),
    tstart_1=tic;
    jVmm = vmm_create(map(alp, train_sets),  ALGS{ii}, params);
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

lamda_arr=[5 20 40 80];
figure(1);
mark_arr=['--go','--b+','-kv','-.r'];
for j=1:length(lamda_arr)
    plot(p_arr,err_arr(:,j), mark_arr(j),'LineWidth',0.8,'MarkerSize',4'); hold on;
end
xlabel('用户离开概率p'); ylabel('预测误差');
legend('\lambda = 5','\lambda = 20','\lambda = 40','\lambda = 80')
 path='D:/doc/PapaerLibrary/Figures/vmm_environment';
 print(path,'-dpng','-r500');

%=======================================================================
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

D_arr=[3 5 10 20 30 40 50];
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

D_arr=[3 5 10 20 30 40 50];
figure(1);
subplot(2,1,1);plot(err_arr,'LineWidth',1.0);
legend('D = 3','D = 5','D = 10','D = 20','D=30','D=40','D=50');
xlabel('训练步数'); ylabel('预测误差');
subplot(2,1,2);plot(train_time_arr,'LineWidth',1.0);
legend('D = 3','D = 5','D = 10','D = 20','D=30','D=40','D=50');
xlabel('训练步数'); ylabel('训练时长/sec');
print('Figs/vmm_oder','-dpng','-r500');


