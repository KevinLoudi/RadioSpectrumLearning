% ****************************************************
% Note: This is a part of the final code aims at spatial locating with RSSI
%  Author: Kevin
%  Date: 15th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 15th January, 2017)
% ****************************************************

%sensor's coordinates
 X=[12.2 15.7 5.2 3.3 11.5 13.4 7.0 4.2];
 %distance to the source
 D=[4.1 5.2 1.3 2.4];
 %real position of the source
 Real1=5.0; Real2=6.3;
 
 cla
circles(X(1), X(5), D(1), 'edgecolor', [0 0 0],'facecolor', 'none','linewidth',4); %AP1 - black
circles(X(2), X(6), D(2), 'edgecolor', [0 1 0],'facecolor', 'none','linewidth',4); %AP2 - green
circles(X(3), X(7), D(3), 'edgecolor', [0 1 1],'facecolor', 'none','linewidth',4); %AP3 - cyan 
circles(X(4), X(8), D(4), 'edgecolor', [1 1 0],'facecolor', 'none','linewidth',4); %AP4 - yellow
axis([0 10 0 10])
hold on
tbl = table(X, d);
d = D.^2;
weights = d.^(-1);
weights = transpose(weights);
beta0 = [5, 5];
modelfun = @(b,X)(abs(b(1)-X(:,1)).^2+abs(b(2)-X(:,2)).^2).^(1/2);
mdl = fitnlm(tbl,modelfun,beta0, 'Weights', weights);
b = mdl.Coefficients{1:2,{'Estimate'}}
scatter(b(1), b(2), 70, [0 0 1], 'filled')
scatter(real1, real2, 70, [1 0 0], 'filled')
hold off