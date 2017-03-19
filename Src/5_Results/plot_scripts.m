% Created on WED Mar 17th 13:25:45 2017
% Propose: Plot scripts 
%
% Enviroment: Matlab 2015b
% @auththor: kevin

%% Plot sumarized spectrum figure

clear; clc; close all;
yes_figure=1;
%FM
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_88_108.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);
freqix=88:0.025:108;
data=dataLevel;
clear dateStamp; clear dataLevel;
display('Data loading finished!!!');

%%  GSM1800
  clear;clc;close all;
 load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_1710_1740.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);
freqix=1710:0.025:1740;
%cut off
f1=(1730-1710)/0.025;
freqix=freqix(f1:end);
data=dataLevel(:,f1:end);
clear dataLevel;
display('GSM 1800 data loaded!!!');

%% TV
clear;clc;close all;
load('D:\\Code\\WorkSpace\\ThesisCode\\Src\\1_Generate_Dataset\\SingalSensorDataset\\Dataset_2_600_798.mat');
timeix=datetime(dateStamp,'InputFormat','yyyy-MM-dd HH:mm:SS');
freqix=600:0.025:798;
t1=datetime('2015-12-16 00:05:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t1_ix=find(timeix==t1);
t2=datetime('2015-12-16 23:55:00','InputFormat','yyyy-MM-dd HH:mm:SS'); t2_ix=find(timeix==t2);

f1=(720-600)/0.025; f2=(750-600)/0.025;
freqix=freqix(:,f1:f2);
data=dataLevel(:,f1:f2);
clear dateStamp; clear dataLevel;
display('TV data loaded!!!');


 %% process missing data
 start_ti=timeix(1);
 step_ti=timeix(2)-timeix(1);
 stop_ti=timeix(end);
 
 ti_slots=(stop_ti-start_ti)/step_ti+1;
 pro_data=0;
 if ti_slots==length(data(:,1))
     display('no missed data');
     pro_data=data;
 %unequal to the therotical value, must miss some data
 else
     display('need to process missed data!!!');
     pro_data=zeros(ti_slots,length(freqix)); %prepare to process the missing data
     ti=start_ti; ixt=1; %index for time slots
     ix=1; %index for the row of dataset
     while ti<= stop_ti
         if  ti~=timeix(ix) %need special treatment
            pro_data(ixt,:)=NaN;
            ixt=ixt+1;
         else %just move the exist data to another matrix
            pro_data(ixt,:)=data(ix,:); 
            ixt=ixt+1; ix=ix+1;
         end
         ti=ti+step_ti;
     end
 end

 %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h_1=figure(1); set(gca,'Position',[.05 .05 .9 .9]);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (FM 88-108MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HHPM')
  %print('Figs/fm','-dpng','-r500');
  display('FM analysised...');
  
  %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h2=figure(3);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (GSM1800 1730-1740MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  axis([1730 1740 0 1.0]);
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HHPM')
  %print('Figs/gsm','-dpng','-r500');
  display('GSM analysised...');
  
  %% do plot
  n_data=Normalize_matrix(pro_data(t1_ix:t2_ix,:),0,1);
  h2=figure(3);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (TV 600-800MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  axis([720 750 0 1.0]);
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HH')
  print('Figs/tv','-dpng','-r500');
  display('TV analysised...');
  
  %% ROHT-ED
  n_data(isnan(n_data(:,1)),:)=[]; %elimate NaN elements
  [cut_point]=Recursive_oneside_hypthesis_testing(n_data, 100);
  thres=cut_point(end);
  %thresholding thresholding
 freqwidth=0.25; %1MHz for FM signal
 [s_cs,s_thres]=SlidingThresholding(n_data,freqix,freqwidth);
 figure(3);  
 %direct threholding
 cs=n_data>thres;
 s_cs=cs|s_cs;
 subplot(3,1,1); imagesc(freqix,time,n_data); xlabel('频率/MHz'); ylabel('时间'); title('频谱监测数据'); datetick('y','HHPM');
 subplot(3,1,2); imagesc(freqix,time,cs); xlabel('频率/MHz'); ylabel('时间'); title('ROHT-能量检测结果'); datetick('y','HHPM');
%  subplot(4,1,3);  plot(s_thres); %title('Adaptive thresholding');%plot(1:length(cut_point),cut_point); title('Thresholding spectrum');
%  title('Sliding-ROHT的阈值');xlabel('窗口（窗宽250kHz）'); ylabel('阈值（已归一化）'); axis([0 81 0 1.5]);
 subplot(3,1,3); imagesc(freqix,time,s_cs); xlabel('频率/MHz'); ylabel('时间'); title('Sliding-ROHT-能量检测结果（窗宽250kHz）'); datetick('y','HHPM');
 %print('Figs/Adoptive thresholding','-dpng');
  
 %% 
 data=Normalize_matrix(pro_data,0,1)>thres;
 time=datenum(timeix(1):step_ti:timeix(end));
 figure(1);
 subplot(3,1,1); imagesc(data); xlabel('频率/MHz'); ylabel('采样次数'); title('（a）频谱占用状态'); %datetick('y','mm/dd HH');
 scr=sum(data,2)/length(freqix); %average each row
 subplot(3,1,2); plot(time(1:1000),scr(1:1000)); title('（b）业务聚集度'); datetick('x',' mm/dd HH','keepticks'); xlabel('月/天 时'); ylabel('聚集度值');
 
 subplot(3,1,3);
 cvd=zeros(size(freqix));
 for freq=1:length(freqix)
    duration=FindZerosBlock(data(:,freq)');
    cvd(freq)=sum(duration)/length(duration);
 end
 cvd_ix=find(cvd>300);
 cvd(cvd_ix)=[];
 histogram(cvd,200,'Normalization','probability'); xlabel('空闲窗口数'); ylabel('频率'); title('（c）持续空闲时长分布');
 %print('Figs/index','-dpng','-r500');
 
 %*********************************************************************************
 %% Test compression ability
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
%% split original data into train-set and test-set
% seq_str='11111111111111111001111111111111111111111111111111';
alp = alphabet(seq_str);
split_point=length(seq_str)-10;
seq_train=seq_str(1:split_point);
seq_test=seq_str((split_point+1):end); 
%algorithm
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
 
disp('---------------------------------------------------');
disp('working with AB={0,1}');
disp('---------------------------------------------------');
disp(' ');

% 3. run each of the VMM algorithms
for i=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{i} ));
    disp('--------')
    %create a VMM through training 
    jVmm = vmm_create(map(alp, seq_train),  ALGS{i}, params);
    disp(sprintf('Pr(0 | 00) = %f', vmm_getPr(jVmm, map(alp,'0'), map(alp,'00'))));
    disp(sprintf('Pr(1 | 00) = %f', vmm_getPr(jVmm, map(alp,'1'), map(alp,'00'))));
    % calculates the length in bits of the  "compressed" representation of
    % seq.  -log[ Pr ( seq | jVmm) ]
    disp(sprintf('-lg(Pr(tar))=%f', vmm_logEval(jVmm,map(alp, seq_test))));
    disp('--------')
    disp(' ');
end
%% % clear; clc; close all;
% len=1e5;
% [s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len);
% seq_str=strjoin(rcs,'');
% seq_mean=sum(rcs)/length(rcs);

%% Prediction ability
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

rand_sel=floor((length(test_sets)-10)*rand());
context=test_sets(rand_sel-10:rand_sel);
rcs=test_sets((rand_sel+1):(rand_sel+10));
pre_steps=length(rcs);
pre_states=zeros(1,pre_steps);


%% 
display('Test prediction ability....');

for ii=1:length(ALGS),
    disp(sprintf('Working with %s', ALGS{ii} ));
    disp('--------')
  for i=1:pre_steps
    jVmm = vmm_create(map(alp, train_sets),  ALGS{ii}, params);
    p_0=vmm_getPr(jVmm, map(alp,'0'), map(alp,context));
    p_1=vmm_getPr(jVmm, map(alp,'1'), map(alp,context));
%     context,
    if p_0>p_1
        pre_states(i)=0;
%         if rcs(i)==num2str(pre_states(i)) display('Correct!!'); end
    else
        pre_states(i)=1;
%         if rcs(i)==num2str(pre_states(i)) display('Correct!!'); end
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

 
  %% Plot simulation data 
t = 0;%time for spectrum access event happen
n = 0;%access times
a = [];%access- time table 
len=1e3;
traffic_data =zeros(1,len);
while(t<len)
    e = exprnd(80);%interval of spectrum access--following expontiental distribution
    t = t + e;
    n = n + 1;
    a = [a;n,t];
end    
a = a'; 

occupy = geornd(0.2,1,n);%following geo-distribution
%set average occupancy time 
for j = 1:n
    i = ceil(a(2,j));
    m = occupy(j);
    traffic_data(i:i+m) = 1;
end 


%make the start point number equal with  stop point
if traffic_data(end)>eps 
    traffic_data=[traffic_data 0.0];
end
traffic=traffic_data>eps;
start_point=find(diff(traffic)==1)+1;
stop_point=find((diff(traffic)==-1));

%randomly pick signal distributions
s_mean_mu=40; s_mean_theata=10;
s_var_mu=10; s_var_theata=3;
s_num=4;
noise_mu=10; noise_theate=3;

mu=random('Normal', s_mean_mu, s_mean_theata, 1, s_num);
theata=random('Normal', s_var_mu, s_var_theata, 1, s_num);
s_ix=floor(s_num*rand(length(start_point),1)+1);
ns=random('Normal',noise_mu,noise_theate,1,length(traffic_data));


figure(1); set(gca,'Position',[.05 .05 .9 .9]); subplot(2,1,1);stem(traffic_data); xlabel('时间'); ylabel('占用状态'); title('（a）频谱接入状态');
axis([0 len eps 1.2]);
hold on;
for i=1:length(s_ix)
    text(start_point(i),1.1,sprintf('PU-%s',int2str(s_ix(i))));
end

traffic_data=ns;

for i=1:length(start_point)
     start_ix=start_point(i);  stop_ix=stop_point(i);
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)+random...
         ('Normal',mu(s_ix(i)),theata(s_ix(i)),1,stop_ix-start_ix+1);  %access from a randomly picked first user
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)-ns(start_ix:stop_ix); %remove the noise part
end

subplot(2,1,2);plot(traffic_data); xlabel('时间'); ylabel('频谱能量'); title('（b）频谱能量变化');
axis([0 len -2*noise_mu 2*s_mean_mu]);

print('Figs/sim','-dpng','-r500');
%% ROHT Process
rcs=traffic;
s=traffic_data;
split=floor(length(s)/2);
train_s=s(1:split);
[cut_point]=Recursive_oneside_hypthesis_testing(train_s, 100);
test_s=s(split+1:end);
test_rcs=rcs(split+1:end);
thres=cut_point(end);
test_ycs=test_s>thres;
err=test_ycs-test_rcs;
      
err_fa=(err==1);
p_fa=sum(err_fa)/length(test_rcs);
      
err_d=(err==0);
p_d=sum(err_d)/length(test_rcs);



%% ROHT-ED with sim
clear; clc; close all;
len=1e3;
exp_time=1e3;
p_fa=zeros(1,exp_time); p_md=zeros(1,exp_time); %detection results of each time
rr=zeros(1,exp_time); %signal-noise ratio
for ii=1:exp_time
    [s,rcs,rr(1,ii)]=Generate_simulation_dataset_v2(42,3,len); %import noise parameters with total length
    split=floor(length(s)/2);
    train_s=s(1:split);
    [cut_point]=Recursive_oneside_hypthesis_testing(train_s, 100);
    
    test_s=s(split+1:end);
    test_rcs=rcs(split+1:end);
    thres=cut_point(end);
    test_ycs=test_s>thres;
    err=test_ycs-test_rcs;
      
    err_fa=(err==1);
    p_fa(ii)=sum(err_fa)/length(test_rcs);
      
    err_d=(err==0);
    p_d(ii)=sum(err_d)/length(test_rcs);
      
    if p_d(ii)<=0.6
          display('abnormal!!!');
    end
      
      err_all=sum(abs(test_ycs-test_rcs))/length(test_rcs);
end
 
 %sum(rr)/exp_time
 display('False Alarm:');
 sum(p_fa)/length(p_fa)
 min(p_fa),max(p_fa)
 display('Detected:');
 sum(p_d)/length(p_d)
 max(p_d),min(p_d)
 figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
 subplot(2,1,2); plot(p_d); xlabel('Times'); ylabel('Probability');title('Detection probability');

%% Plot energy-detection ROC curve
clear; clc; close all;
snr=[-8,-9,-10,-11];
num=length(snr);
figure(1);
ix=['（a）','（b）','（c）','（d）'];
color={'--g','--b','--k','--r'};
color2={'g','b','k','r'};

figure(1)
for i=1:length(snr)
    [Pf,Pd,Pd_the]=Energy_detection_simulation(snr(i));
    subplot(2,2,i);
    plot(Pf, Pd,color{i},'LineWidth',2,'MarkerEdgeColor',color2{i}, 'MarkerSize',6); hold on; %paractical
    plot(Pf, Pd_the, color2{i},'LineWidth',2); %theory
    title(sprintf('SNR=%s dB',snr(i)));
    xlabel({['P_{fa}'],[ix(i)]}); ylabel('P_{d}');
    legend('实验值','理论值');
end
xlabel('P_{fa}'); ylabel('P_{d}');
%legend('SNR=-8dB','SNR=-9dB','SNR=-10dB','SNR=-11dB');
axis([0 0.5 0.5 1]);

print('Figs/ed','-dpng','-r500');

%% 

  
  
