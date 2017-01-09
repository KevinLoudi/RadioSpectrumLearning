clear;
load level_day1.mat; %dlp 80-108MHz
dlp_m=uint8(normalize_matrix(dlp,0,255));
%imagesc(dlp_m);
imshow(dlp_m);
imwrite(dlp_m,'level_day1.jpg');