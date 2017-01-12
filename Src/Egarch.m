% ****************************************************
% Note: This is a part of the final code aims at ARIMA time-series fitting
%   with exponential generalized autoregressive conditional 
%   heteroscedastic (EGARCH) model
%  Author: Kevin
%  Date: 12th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 12th January, 2017)
% ****************************************************

%%EGARCH conditional variance time series model
clear; clc;
%create  EGARCH(p,q) model
p=1; q=1;
% conditional mean model offset is zero
Mdl = egarch('GARCHLags',p,'ARCHLags',q,'LeverageLags',1,'Offset',NaN);

%load data and plot
load Data_Danish;
nr = DataTable.RN;
figure(1);
plot(dates,nr);
hold on;
plot([dates(1) dates(end)],[0 0],'r:'); % Plot y = 0
hold off;
title('Danish Nominal Stock Returns');
ylabel('Nominal return (%)');
xlabel('Year');

%estimate model with the input time-series
EstMdl = estimate(Mdl,nr);

%simulate 
numObs = numel(nr); % Sample size (T)
numPaths = 100;     % Number of paths to simulate
rng(1);             % For reproducibility
[VSim,YSim] = simulate(EstMdl,numObs,'NumPaths',numPaths);

%make comparsion
VSimBar = mean(VSim,2);
VSimCI = quantile(VSim,[0.025 0.975],2);
YSimBar = mean(YSim,2);
YSimCI = quantile(YSim,[0.025 0.975],2);

figure;
subplot(2,1,1);
h1 = plot(dates,VSim,'Color',0.8*ones(1,3));
hold on;
h2 = plot(dates,VSimBar,'k--','LineWidth',2);
h3 = plot(dates,VSimCI,'r--','LineWidth',2);
hold off;
title('Simulated Conditional Variances');
ylabel('Cond. var.');
xlabel('Year');

subplot(2,1,2);
h1 = plot(dates,YSim,'Color',0.8*ones(1,3));
hold on;
h2 = plot(dates,YSimBar,'k--','LineWidth',2);
h3 = plot(dates,YSimCI,'r--','LineWidth',2);
hold off;
title('Simulated Nominal Returns');
ylabel('Nominal return (%)');
xlabel('Year');
legend([h1(1) h2 h3(1)],{'Simulated path' 'Mean' 'Confidence bounds'},...
    'FontSize',7,'Location','NorthWest');

%make forecast for future 10 steps
numPeriods = 10;
vF = forecast(EstMdl,numPeriods,'Y0',nr);
v = infer(EstMdl,nr);

figure;
plot(dates,v,'k:','LineWidth',2);
hold on;
plot(dates(end):dates(end) + 10,[v(end);vF],'r','LineWidth',2);
title('Forecasted Conditional Variances of Nominal Returns');
ylabel('Conditional variances');
xlabel('Year');
legend({'Estimation sample cond. var.','Forecasted cond. var.'},...
    'Location','Best');
