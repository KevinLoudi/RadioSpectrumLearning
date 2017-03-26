% Created on Fri 17th Mar  13:25:45 2017
% Propose: Generate simulation data according to the 3.4
%   detection and thresholding ability
% Enviroment: Matlab 2015b
% @auththor: kevin

function [traffic_data,traffic,r]=Generate_simulation_dataset_v2(noise_mu,noise_theate,len,lamda,p,s_num)
t = 0;%time for spectrum access event happen
n = 0;%access times
a = [];%access- time table 
traffic_data =zeros(1,len);
while(t<len)
    e = exprnd(lamda);%interval of spectrum access--following expontiental distribution
    t = t + e;
    n = n + 1;
    a = [a;n,t];
end    
a = a'; 

occupy = geornd(p,1,n);%following geo-distribution
%set average occupancy time 
for j = 1:n
    i = ceil(a(2,j));
    m = occupy(j);
    traffic_data(i:i+m) = 1;
end 
%% 
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
% s_num=4;

mu=random('Normal', s_mean_mu, s_mean_theata, 1, s_num);
theata=random('Normal', s_var_mu, s_var_theata, 1, s_num);
s_ix=floor(s_num*rand(length(start_point),1)+1);
ns=random('Normal',noise_mu,noise_theate,1,length(traffic_data));
traffic_data=traffic_data+ns;

for i=1:length(start_point)
%      s_mu=random('normal',singal_mu,10);
%      s_theta=random('normal',singal_theate,3);
     start_ix=start_point(i);  stop_ix=stop_point(i);
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)+random...
         ('Normal',mu(s_ix(i)),theata(s_ix(i)),1,stop_ix-start_ix+1);  %access from a randomly picked first user
     traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix)-ns(start_ix:stop_ix); %remove the noise part
%      traffic_data(start_ix:stop_ix)=traffic_data(start_ix:stop_ix).*...
%      random('normal',s_mu,s_theta,1,stop_ix-start_ix+1);
end

r=log10(s_mean_mu/noise_mu);
%r=snr(traffic_data-ns,ns);
%ix=find(traffic_data<eps); traffic_data(ix)=eps;


end