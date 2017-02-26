% Created on Mon 26th Feb  13:25:45 2017
% Propose: Thresholding for the input data
% Enviroment: Matlab 2015b
% @auththor: kevin
% ThresInfo.LevelArr, Level, Mask

function [ThresInfo]=DoubleThresholding(Data)
    maxIter=100; tol=1e-6;
    img=Data;
    delete Data;
    imgFlat = img(:);
    ii = 1; %iteration index
    thresh(ii) = mean(imgFlat); %threshold of the ii-th iteration
    %start iteration
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
    
    %save compuation results
    ThresInfo.LevelArr=thresh;
    ThresInfo.Level=thresh(end);
    ThresInfo.Mask=imgBw;
end