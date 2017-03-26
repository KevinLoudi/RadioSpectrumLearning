 clear; clc; close all;
 s_mask=zeros(1e3,1e3);
 %define some channel
 s_mask(: ,[50:90, 140:185, 214:238, 290:310, 460:480, 665:675, 730:740, 800:835, 890:924])...
     =s_mask(: ,[50:90, 140:185, 214:238, 290:310, 460:480, 665:675, 730:740, 800:835, 890:924])|1;
 
 c_ix={50:90, 140:185, 214:238, 290:310, 460:480, 665:675, 730:740, 800:835, 890:924};
 c_num=9;
 
 spectrum=s_mask;
 len=1e3;
 
 %generate mimic spectrum
 figure(1);
 for j=1:c_num
         % s_num:4  lamda:80 p=0.5
        [s,rcs,rr]=Generate_simulation_dataset_v2(10,3,len,100,0.2,4); %import noise parameters with total length
    for i=c_ix{j}
        spectrum(:,i)=spectrum(:,i).*s(1,1:1e3)';
    end
 end
 imagesc(spectrum); xlabel('频率'); ylabel('时间');title('轻负载');
 colorbar;
 
path='D:/doc/PapaerLibrary/Figures/data_mul';
print(path,'-dpng','-r500');
 


 
% freqix=1:1e3; freqwidth=100;
% timeix=1:1e3; timewidth=100;
% [s_cs,s_thres]=SlidingThresholding(spectrum,freqix,freqwidth);

