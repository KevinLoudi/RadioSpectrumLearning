% Created on Mon 26th Feb  13:25:45 2017
% Propose: Hist analysis after thresholding 
% Enviroment: Matlab 2015b
% @auththor: kevin
% ThresInfo.LevelArr, Level

function [HistInfo]=HistAnalysis(Data, Threshold)
    figure;
    c=Data(:); %reshape input data
    c1=c(c<Threshold); %select sample data that belongs to specific class
    c2=c(c>Threshold);
    h1=histogram(c1,'Normalization','probability');
    pd1 = fitdist(c1,'Normal'); %fit sample to normal distribution
    hold on;
    h2=histogram(c2,'Normalization','probability');
    pd2 = fitdist(c2,'Normal');
    
    %plot probablity distribution
    x1.min=min(c); x1.max=Threshold;
    x1.value=x1.min:((x1.max-x1.min)/100):x1.max
    y1.value=pdf(pd1,x1.value);
    
    x2.min=Threshold; x2.max=max(c);
    x2.value=x2.min:((x2.max-x2.min)/100):x2.max
    y2.value=pdf(pd2,x2.value);
    hold on;
    plot(x1.value, y1.value, 'LineWidth',2);
    hold on;
    plot(x2.value, y2.value, 'LineWidth',2);
    title('classify analysis');
    
    HistInfo.LT=pd1.mu+2*pd1.sigma;
    HistInfo.HT=pd2.mu-2*pd2.sigma;
end
