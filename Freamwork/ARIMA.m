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
AR=[1 0.5];  %AR coeffs
MA=[1 0.4 0.3]; %MA coeffs
T=1000;
WN=randn(T,1); %generate gaussian white noise
y=filter(MA,AR,WN); %generate y

fvtool(MA,AR,'impulse'); %create impulse response

[c lags]=xcov(y,'biased');
figure;
plot(lags,c);
title('Sample Convariance');

