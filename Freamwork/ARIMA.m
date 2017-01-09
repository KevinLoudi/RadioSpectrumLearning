%*******************************************************
% Author: Kevin
% Propose: ARIMA model Matlab version
% Date: 29th Dec, 2016
% Environment: Matlab 2015b
%*******************************************************

% clear; clc; clf;
% %specific a nonseasonal ARIMA model
% %two nonseasonal AR coefficients ( $p$ = 2), 
% %two nonseasonal MA coefficients ( $q$ = 2), 
% %and one degree of differencing ( $D$ = 1)
% Mdl = arima(2,1,2);
% 
% %Load the NASDAQ data included with the toolbox. 
% %Extract the first 1500 observations of the Composite Index
% load Data_EquityIdx
% nasdaq = DataTable.NASDAQ(1:1500);
% 
% %Fit the model to the first half of the data.
% EstMdl = estimate(Mdl,nasdaq(1:750));

clear;
%**************Simulating an ARMA Model
AR=[1 0.5];  %AR coeffs
MA=[1 0.4 0.3]; %MA coeffs
T=1000;
WN=randn(T,1); %generate gaussian white noise
%  a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na) filter(B,A,X)
y=filter(MA,AR,WN); %generate data from ARMA model

%graph the impulse response of an ARMA
% computes the
%     Magnitude Response for the filter defined by numerator and denominator
%     coefficients in vectors B and A
fvtool(MA,AR,'impulse'); %create impulse response

%Sample Covariances
[c lags]=xcov(y,'biased');
figure;
plot(lags,c);
title('Sample Convariance');

%make into density
w=0:0.1:pi;
h=freqz(MA,AR,w);
sd=abs(h).^2./sqrt(2*pi);

%estimate AR(8)
[sdc wc]=pcov(y,8);
[sdy wy]=pyulear(y,8); %using the Yule?walker equations
[sdp wp] = periodogram(y,[],'onesided'); % estimate using sample
