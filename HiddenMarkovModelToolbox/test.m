%*******************************************************
% Author: Kevin
% Propose: Testing for HMM Toolbox
% Date: 29th Dec, 2016
% Environment: Matlab 2015b
%*******************************************************

% %test for all build-in method
% clear; clc; clf;
% d = 2; %internal state numbers
% k = 2; %visible 
% n = 100; %length
% [x, model] = hmmRnd(d, k, n);
% %estimate internal states through observances and model
% z = hmmViterbi(x,model);
% 
% figure(1);
% plot(z,'b');
% %calculate the probability of occuring the specific seq
% %with logically likehood
% %forward filtering algorithm
% [alpha,llh] = hmmFilter(x,model);
% 
% %EM algorithm to fit the parameters of HMM model 
% [model_est, llh] = hmmEm(x,k);
% z_est=hmmViterbi(x,model_est);
% figure(2)
% % plot(z,'b'); hold on;
% % plot(z_est,'r');
% plot(z-z_est,'r');
% err=sum(z-z_est)/length(z),
% %plot(llh)

%-------------------------------------------------------------
clear;clc;clf;
%simulate a hidden state seq
A = [0.95,0.05;0.10,0.90]; %two internal state
B = [1/6, 1/2 2/6;1/2, 1/10, 2/5;]; %three visible state
seq = hmmgenerate(1000,A,B);

%estimate transfer matrix
[AHat, BHat] = hmmtrain(seq,A,B);
%state prodiction
pStates = hmmdecode(seq,A,B);
pStates_est=hmmdecode(seq,AHat, BHat);
figure(1)
plot(pStates(1,:),'b'); 
hold on;
plot(pStates_est(1,:),'g');



