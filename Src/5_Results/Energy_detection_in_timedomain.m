%Energy Detection in the Time Domain

FrameLength = 20;
Fs = 100;
movrmsWin = dsp.MovingRMS(20);
scope  = dsp.TimeScope('SampleRate',Fs,...
    'TimeSpanOverrunAction','Scroll',...
    'TimeSpan',100,...
    'ShowGrid',true,...
    'YLimits',[-1.0 350],'LayoutDimensions',[3 1],'NumInputPorts',3);

scope.ActiveDisplay = 1;
scope.YLimits = [0 5];
scope.Title = 'Input Signal';

scope.ActiveDisplay = 2;
scope.Title = 'Compare Signal Energy with a Threshold';

scope.ActiveDisplay = 3;
scope.YLimits = [0 2];
scope.PlotType = 'Stairs';
scope.Title = 'Detect When Signal Energy Is Greater Than the Threshold';

count = 1;
Vect = [1/8 1/2 1 2 3 4 3 2 1];
index = 1;
threshold = 200;
for index = 1:length(Vect)
    V = Vect(index);
    for i = 1:80
        x = V + 0.1 * randn(FrameLength,1);
        y1 = movrmsWin(x);
        y1ener = (y1(end)^2)*20;
        event = (y1ener>threshold);
        scope(y1,[y1ener,threshold],event);
    end
end