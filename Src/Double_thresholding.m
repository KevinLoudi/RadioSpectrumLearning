%% load data and do iterative threshold
clear; clc;
[x,map]=imread('cameraman.tif');

img=x;
maxIter = 1e2;
tol = 1e-6;
imgFlat = img(:);
ii = 1; %iteration index
thresh(ii) = mean(imgFlat); %threshold of the ii-th iteration
while ii<maxIter
    display('=======================');
    display('New iteration....');
	imgBw = img>thresh(ii);
	mbt = mean(imgFlat(~imgBw)); %mean of the backgroud class
    vbt = std(double(imgFlat(~imgBw)));
    
	mat = mean(imgFlat(imgBw)); %mean of the foreground class
    vat = std(double(imgFlat(imgBw)));
    
    if thresh(ii)>mbt+3*vbt
        display('fine threshold for background');
    elseif  thresh(ii)<mat-3*vat   
        display('fine threshold for foreground');
    end
    
	thresh(ii+1) = 0.5*(mbt + mat); %put the threshold inside two mean
	if abs(thresh(ii+1) - thresh(ii)) > tol
		ii = ii + 1;
	else
		break
	end;
end;

level = thresh(end);
%% plot histogram of the two classified group and fit into normal distributions

clf;
%plot sample of two classes separtely
s1=imgFlat(~imgBw);
h1=histogram(s1,'Normalization','probability');
hold on;

%fit distribution and plot pdf
pd1 = fitdist(s1,'Normal');
x1.min = 0;%pd1.mu-3.*(pd1.sigma);
x1.max = level;%pd1.mu+3*pd1.sigma;
%plot fitted distribution
x1.value=x1.min:((x1.max-x1.min)/100):x1.max
y1.value=pdf(pd1,x1.value);
hold on;
plot(x1.value, y1.value, 'LineWidth',2);

%plot sample of two classes separtely
s2=imgFlat(imgBw);
h2=histogram(s2,'Normalization','probability');

%fit distribution and plot pdf
pd2 = fitdist(s2,'Normal');
x2.min = level;%pd2.mu-3.*(pd2.sigma);
x2.max = 255;%pd2.mu+3*pd2.sigma;
%plot fitted distribution
x2.value=x2.min:((x2.max-x2.min)/100):x2.max
y2.value=pdf(pd2,x2.value);
hold on;
plot(x2.value, y2.value, 'LineWidth',2);


%% double-threshold, separate into three region
thres_1=x1.max;
thres_2=x2.max

