
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