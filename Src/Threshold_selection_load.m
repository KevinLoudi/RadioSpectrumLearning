% ****************************************************
%  Note: This is a part of the final code aim at select a threshold for
%      energy detection, return threshold vector
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************
function thresh=Threshold_selection_load(r)
    %load data
    level=r;

    %Global image threshold selection
    thresh = multithresh(level,10);
    %imwrite(level,'../Data/level.jpg');
end