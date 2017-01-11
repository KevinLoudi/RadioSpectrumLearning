%read image
he = imread('level_day1.jpg');
he = ind2rgb(he,colormap);
imshow(he), title('Spectrum image');
text(size(he,2),size(he,1)+15,'Spectrum of Kevin, Sichuan University', ...
     'FontSize',7,'HorizontalAlignment','right');
 
 %convert image from RGB to L*a*b*color
cform = makecform('srgb2lab');
lab_he = applycform(he,cform);
imshow(lab_he);

%classify colors in a*b space 
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);

nColors = 2;
nClasses = 2;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',nClasses);
                                  
%label every pixel                                  
segmented_images = cell(1,nClasses);
rgb_label = repmat(pixel_labels,[1 1 nClasses]);

for k = 1:nColors
    color = he;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
end

%show the segemented image
figure(1);
imshow(segmented_images{1}), title('objects in cluster 1');
figure(2);
imshow(segmented_images{2}), title('objects in cluster 2');
% figure(3);
% imshow(segmented_images{3}), title('objects in cluster 3');

