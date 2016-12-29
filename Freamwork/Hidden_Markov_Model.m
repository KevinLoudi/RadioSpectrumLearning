%*******************************************************
% Author: Kevin
% Propose: Implement HMM in channel state modeling
% Date: 28th Dec, 2016
% Environment: Matlab 2015b
%*******************************************************

clf;
clear;
clc;
%-----------test part
%-----------Generating a Test Sequence
TRANS = [.9 .1; .05 .95;]; %transition matrix 
EMIS = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6; 7/12, 1/12, 1/12, 1/12, 1/12, 1/12]; %emissions matrix
%generate a random sequence of states and emissions from the model
%output seq is the sequence; output states is the sequence of states
[seq,states] = hmmgenerate(1000,TRANS,EMIS);
figure(1),
plot(seq,'g');
hold on;
plot(states,'r');

%-------------Estimating the State Sequence
%uses the Viterbi algorithm to compute the most likely sequence of states
%likelystates is a sequence the same length as seq
likelystates = hmmviterbi(seq, TRANS, EMIS);
%test the accuracy of hmmviterbi
accuracy=sum(states==likelystates)/1000;
accuracy,
figure(2),
% plot(states,'r');
% hold on;
% plot(likelystates,'b');
plot(states-likelystates,'b');

%-------------Estimating Transition and Emission Matrices
%estimate the transition and emission matrices TRANS and EMIS 
%given a sequence seq of emissions
tic;
[TRANS_EST, EMIS_EST] = hmmestimate(seq, states); %Using hmmestimate

%have initial guesses for TRANS and EMIS, you can still estimate TRANS and
%EMIS using hmmtrain
TRANS_GUESS = [.85 .15; .1 .9];
EMIS_GUESS = [.17 .16 .17 .16 .17 .17;.6 .08 .08 .08 .08 08];
%hmmtrain uses an iterative algorithm that alters the matrices TRANS_GUESS 
%and EMIS_GUESS so that at each step the adjusted matrices are more likely to 
%generate the observed sequence, seq
tol=10;
[TRANS_EST2, EMIS_EST2] = hmmtrain(seq, TRANS_GUESS, EMIS_GUESS,'tolerance', tol)
running_time=toc;
%two factors reduce the reliability of hmmtrain
%converge to a local maximum 

%-------------Estimating Posterior State Probabilities
% the conditional probabilities that the model is in a particular state when 
% it generates a symbol in seq
%The output PSTATES is an M-by-L matrix, where M is the number of states 
%and L is the length of seq. PSTATES(i,j) is the conditional probability that the 
%model is in state i when it generates the jth symbol of seq
[PSTATES,logpseq] = hmmdecode(seq,TRANS,EMIS);%logarithm of the probability of the sequence seq

%The probability of a sequence tends to 0 as the length of the sequence increases
