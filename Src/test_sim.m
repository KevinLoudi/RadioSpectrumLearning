   clear; clc; close all;
   orign_path='D:/Code/WorkSpace/ThesisCode/Src/4_Spatio_time/SpatialDataset/%s';
   load(sprintf(orign_path, 'Spdata_1068_1068.mat'));
   load(sprintf(orign_path, 'Spdata_1068_1068.mat'));
   load(sprintf(orign_path,'Splocation_1068_1068.mat'));
   % create random field with autocorrelation
   %[X,Y] = meshgrid(locateion_sp(1,1:end),locateion_sp(1:end,1));
   Z=sum(data_sp,1,'omitnan')./1e6;
   v = variogram([X(1,:) Y(:,1)],Z,'plotit',false,'maxdist',100);
    %% 
    
    Z = randn(size(X));
    Z = imfilter(Z,fspecial('gaussian',[40 40],8));

    %% 
    
    % sample the field
    n = 500;
    x = rand(n,1)*500;
    y = rand(n,1)*500;
    z = interp2(X,Y,Z,x,y);
    
      % plot the random field
    subplot(2,2,1)
    imagesc(X(1,:),Y(:,1),Z); axis image; axis xy
    hold on
    plot(x,y,'.k')
    title('random field with sampling locations')

    % calculate the sample variogram
    v = variogram([x y],z,'plotit',false,'maxdist',100);
    % and fit a spherical variogram
    subplot(2,2,2)
    [dum,dum,dum,vstruct] = variogramfit(v.distance,v.val,[],[],[],'model','stable');
    title('variogram')

    % now use the sampled locations in a kriging
    [Zhat,Zvar] = kriging(vstruct,x,y,z,X,Y);
    subplot(2,2,3)
    imagesc(X(1,:),Y(:,1),Zhat); axis image; axis xy
    title('kriging predictions')
    subplot(2,2,4)
    contour(X,Y,Zvar); axis image
    title('kriging variance')






 