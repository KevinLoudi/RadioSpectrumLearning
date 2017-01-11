% ****************************************************
% Note: This is a part of the final code aim at deciding specific channel
%    state
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

function channe_status=Decide_status()

%get threshold vector
thresh=Threshold_selection();

%select a threshold
a_thresh=thresh(5); %threshold=98

%decide channel status
channe_status=level>a_thresh;

%plot channe status
imagesc(channe_status);

end





