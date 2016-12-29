%*******************************************************
% Author: Kevin
% Propose: Testing for HMM Toolbox
% Date: 29th Dec, 2016
% Environment: Matlab 2015b
%*******************************************************

%test for all build-in method
clear; clc; clf;
d = 2; %internal state numbers
k = 2; %visible 
n = 100; %length
[x, model] = hmmRnd(d, k, n);
%estimate internal states through observances and model
z = hmmViterbi(x,model);

figure(1);
plot(z,'b');
%calculate the probability of occuring the specific seq
%with logically likehood
%forward filtering algorithm
[alpha,llh] = hmmFilter(x,model);

%EM algorithm to fit the parameters of HMM model 
[model_est, llh] = hmmEm(x,k);
plot(llh)


