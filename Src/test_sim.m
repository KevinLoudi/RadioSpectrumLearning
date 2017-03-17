
clear; clc; close all;
len=1e3;
exp_time=1e3;
p_fa=zeros(1,exp_time); p_md=zeros(1,exp_time); %detection results of each time
for ii=1:exp_time
    [s,rcs]=Generate_simulation_dataset_v2(30,15,5,10,len);
    [cut_point]=Recursive_oneside_hypthesis_testing(s, 100);
    
     thres=cut_point(end);
      ycs=s>thres;
      err=ycs-rcs;
      
      err_fa=(err==1);
      p_fa(ii)=sum(err_fa)/length(rcs);
      
      err_md=(err==-1);
      p_md(ii)=sum(err_md)/length(rcs);
      
      err_all=sum(abs(ycs-rcs))/length(rcs);
end

 sum(p_fa)/length(p_fa)
 sum(p_md)/length(p_md)
 figure(7); subplot(2,1,1); plot(p_fa); xlabel('Times'); ylabel('Probability');title('False alarm probability');
 subplot(2,1,2); plot(p_md); xlabel('Times'); ylabel('Probability');title('Miss detection probability');


