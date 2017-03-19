%% 
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
  h2=figure(3);
  subplot(2,1,1); plot(freqix,n_data(1,:)); title('2015-12-16 00:05:00 (GSM1800 1730-1740MHz)'); ylabel('频谱能量(归一化)/dB\muV^{-1}'); xlabel('频率/MHz');
  axis([1730 1740 0 1.0]);
  time=datenum(timeix(t1_ix:t2_ix));
  subplot(2,1,2); imagesc(freqix,time,n_data); title('2015-12-16'); ylabel('时间'); xlabel('频率/MHz');
  datetick('y','HHPM')
  %print('Figs/gsm','-dpng','-r500');
  display('GSM analysised...');
  
    %% ROHT-ED
  n_data(isnan(n_data(:,1)),:)=[]; %elimate NaN elements
  figure(2); %histogram(n_data(180:280,:),50);
  [cut_point]=Recursive_oneside_hypthesis_testing(n_data, 100);
  for i=1:length(cut_point)
      n_data(find(n_data>cut_point(i)))=[];
      pd = fitdist(n_data(:),'Normal');
      x=0:0.01:1;
      y=Normalize_matrix(pdf(pd,x),0,1);
      if i==1 | i==3 | i==6 | i==length(cut_point)
         plot(x, y, 'LineWidth',2);
      end
      hold on;
  end
  xlabel('频谱能量（已归一化）');
  legend('迭代1','迭代3','迭代6','迭代完成');
  ylabel('拟合直方图概率');
  print('Figs/ROHT_reduce','-dpng','-r500');
  thres=cut_point(end);
  %thresholding thresholding
 freqwidth=0.25; %1MHz for FM signal
 [s_cs,s_thres]=SlidingThresholding(n_data,freqix,freqwidth);

 
 cs=n_data>thres;
 s_cs=cs|s_cs;
 
 figure(3);
 plot(cut_point); xlabel('迭代次数'); ylabel('阈值');
 print('Figs/roht_sel','-dpng','-r500');
 




