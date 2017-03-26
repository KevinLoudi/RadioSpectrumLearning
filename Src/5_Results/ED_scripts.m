%% ROHT with confiedence level
clear; clc; close all;
confidence_arr=[0.61, 0.69, 0.73, 0.79, 0.82, 0.84, 0.88, 0.90, 0.92, 0.93, 0.95, 0.96,0.97,0.98, 0.986,0.989...
    ,0.992,0.995,0.998,0.999,0.9997];
z_alph_arr=[0.3, 0.5, 0.6, 0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,  1.9, 2.1, 2.2, 2.3...
    ,2.4,2.6,2.9,3.1, 3.4];
len=1e5; confid_len=length(z_alph_arr);
exp_time=1;
p_fa=zeros(confid_len,exp_time); p_md=zeros(confid_len,exp_time); %detection results of each time
rr=zeros(1,exp_time); %signal-noise ratio
for j=1:length(confidence_arr)
    z_alph=z_alph_arr(j);
for ii=1:exp_time
    [s,rcs,rr(1,ii)]=Generate_simulation_dataset_v2(10,3,len,20,0.05,4); %import noise parameters with total length
    split=floor(length(s)/2);
    train_s=s(1:split);
    [cut_point]=Recursive_oneside_hypthesis_testing(train_s, 100,z_alph);
    
    test_s=s(split+1:end);
    test_rcs=rcs(split+1:end);
    thres=cut_point(end);
    test_ycs=test_s>thres;
    err=test_ycs-test_rcs;
      
    err_fa=(err==1);
    p_fa(j,ii)=sum(err_fa)/length(test_rcs);
      
    err_d=(err==-1);
    p_md(j,ii)=sum(err_d)/length(test_rcs);
      
    if p_md(j,ii)<=0.6
          display('abnormal!!!');
    end
      
     err_all=sum(abs(test_ycs-test_rcs))/length(test_rcs);
end
 
end
%% 
 display('False Alarm:');
 sum(p_fa)/length(p_fa)
 min(p_fa),max(p_fa)
 display('Detected:');
 sum(p_md)/length(p_md)
 max(p_md),min(p_md)
 figure(7);  plot(confidence_arr, p_fa,'--rs', 'LineWidth', 1.0); xlabel('置信度', 'FontSize',12); ylabel('概率指标', 'FontSize',12);
 hold on; plot(confidence_arr,p_md,'-b+','LineWidth', 1.0); xlabel('置信度', 'FontSize',12); ylabel('概率指标', 'FontSize',12);
 legend('误警概率','漏检概率', 'FontSize',12); title('重负载，固定占用')
 axis([0.75 1 0 1]);
 path='D:/doc/PapaerLibrary/Figures/ed_confidence';
%  print(path,'-dpng','-r500');
 
 %% ROHT-ED with window size
 clear; clc; close all;
 s_mask=zeros(1000,1000);
 %define some channel
 s_channel=xor(s_mask(: ,50:90, 140:185, 214:238, 290:310, 460:480, 665:675, 730:740, 800:835, 890:924),1);
 
 %% Sim-data
 len=1e3;
  n=4;  lamda=80; p=0.2; n_mu=10; n_theata=3; %light load dc=0.08 snr=0.6
%   n=4; lamda=5; p=0.8; n_mu=10; n_theata=3; %heave load dc=0.14 snr=0.6
%   n=4; lamda=20; p=0.05; n_mu=10; n_theata=3; %heave load, static occupay
%   n=8; lamda=5; p=0.8; n_mu=10; n_theata=3; %heave load, random occupay
%   n=4; lamda=5; p=0.8; n_mu=35; n_theata=3;

 [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 plot(s); hold on; %plot(rcs*100);
 dc=sum(rcs)/len,
 snr,
 
 %% Plot sim-data figures
 
 figure(2)
 subplot(2,2,1) 
 n=4;  lamda=80; p=0.2; n_mu=10; n_theata=3; [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 plot(s); xlabel('时间','FontSize',12); ylabel('频谱能量','FontSize',12); title('轻负载','FontSize',12); axis tight;
 
   subplot(2,2,2) 
 n=4; lamda=20; p=0.05; n_mu=10; n_theata=3; [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 plot( s); xlabel('时间','FontSize',12); ylabel('频谱能量','FontSize',12); title('重负载，固定占用','FontSize',12); axis tight;
 
   subplot(2,2,3) 
 n=8; lamda=5; p=0.8; n_mu=10; n_theata=3; [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 plot(s); xlabel('时间','FontSize',12); ylabel('频谱能量','FontSize',12); title('重负载，随机占用','FontSize',12); axis tight;
 
   subplot(2,2,4) 
 n=4; lamda=5; p=0.8; n_mu=35; n_theata=3; [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 plot(s); xlabel('时间','FontSize',12); ylabel('频谱能量','FontSize',12); title('重负载，低信噪比','FontSize',12);  axis tight; 
 
  path='D:/doc/PapaerLibrary/Figures/four_sim_spectrum';
%  print(path,'-dpng','-r500');

%% OTSU
 clear; clc; close all;
 len=1e2;
 exp_time=1e3;
    n=4;  lamda=80; p=0.2; n_mu=10; n_theata=3; %light load dc=0.08 snr=0.6
%      n=4; lamda=20; p=0.05; n_mu=10; n_theata=3; %heave load, static occupay
%    n=8; lamda=5; p=0.8; n_mu=10; n_theata=3; %heave load, random occupay
%    n=4; lamda=5; p=0.8; n_mu=35; n_theata=3;
 
 for i=1:exp_time

 [s,rcs,snr]=Generate_simulation_dataset_v2(n_mu,n_theata,len,lamda,p,n);
 
 otsu_res=otsu(s,2)-1;
 otsu_err=otsu_res-rcs;
 otsu_err_fa=(otsu_err==1);
 p_fa_otsu(i)=sum(otsu_err_fa)/length(rcs);
 otsu_err_md=(otsu_err==-1);
 p_md_otsu(i)=sum(otsu_err_md)/length(rcs);
 metric_otsu(i)= (1-p_md_otsu(i))/(p_fa_otsu(i)+eps);
%  display('Probability metrics of OTSU: P_fa/P_md '); p_fa_otsu,p_md_otsu, (1-p_md_otsu)/(p_fa_otsu+eps),

 
 z_alph=1.4;
 [cut_point]=Recursive_oneside_hypthesis_testing(s, 100, z_alph);

roht_res=s>cut_point(end);
roht_err=roht_res-rcs;
roht_err_fa=(roht_err==1);
p_fa_roht(i)=sum(roht_err_fa)/length(rcs);
roht_err_md=(roht_err==-1);
p_md_roht(i)=sum(roht_err_md)/length(rcs);
metric_roht(i)=(1-p_md_roht(i))/(p_fa_roht(i)+eps);
% display('Probability metrics of ROHT: P_fa/P_md ');p_fa_roht,p_md_roht, (1-p_md_roht)/(p_fa_roht+eps),

 end
 display('Probability metrics of OTSU: P_fa/P_md ');
 display('Averge:'), sum(p_fa_otsu)/exp_time, sum(p_md_otsu)/exp_time, sum(metric_otsu)/exp_time, 
 display('===================');
 display('Best:'), min(p_fa_otsu), min(p_md_otsu), max(metric_otsu),
 display('===================');
 display('Worst:'),max(p_fa_otsu), max(p_md_otsu), min(metric_otsu),
 display('===================');
 
 display('Probability metrics of ROHT: P_fa/P_md ');
 display('Averge:'), sum(p_fa_roht)/exp_time, sum(p_md_roht)/exp_time, sum(metric_roht)/exp_time, 
 display('===================');
 display('Best:'), min(p_fa_roht), min(p_md_roht), max(metric_roht),
 display('===================');
 display('Worst:'),max(p_fa_roht), max(p_md_roht), min(metric_roht),
 display('===================');
 