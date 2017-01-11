% ****************************************************
% Note: This is a part of the final code aim at loading collected data
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

function level=Load_data()

%load level data and normalize to [0,255]
load ../Data/level_1800; %dlp 80-108MHz
level=uint8(Normalize_matrix(dataLevel(1:200,1:800),0,255));

%pack data into grey image
imshow(level);
imwrite(level,'../Data/level.jpg');

end



