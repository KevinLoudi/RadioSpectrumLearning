% ****************************************************
%  Note: Test script for the project
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

clear;clf;clc;
% %generate signal
% [data,rcs]=Generate_signal();
% %select threshold
% thresh=Threshold_selection_load(data);
% cs=data>thresh(1);
% %error rate 
% rcs=rcs>0;
% err_s=abs(rcs-cs);
% err=sum(abs(rcs-cs))/length(data);
% 
% plot(1:length(data),err_s*100,'b');
% hold on;
% plot(1:length(data),data,'g');
% % hold on;
% % plot(1:length(data),(rcs>0)*100,'k');

 % spatial data represent and analysis in a random field
  % build a struct containing the correlation information
  corr.name = 'gauss';  %correlation type
  corr.c0 = 1; %scaling parameters
  corr.sigma = 1; %variance scaling parameter

  mesh = linspace(-1,1,101)';              % generate a mesh
  data.x = [-1; 1]; data.fx = [0; -1];    % specify boundaries

  % KL: components of the Karhunen-Loeve representation of the random field
  %include mean, eigenvectors covariance matrix and square root of the eigenvalues
  % F: A matrix of random field realizations 
  [F,KL] = randomfield(corr, mesh, ...
              'nsamples', 10, ...
              'data', data, ...
              'filter', 0.95);

  % to generate 100 more samples using the KL
  trunc = length(KL.sv);                  % get the truncation level
  W = randn(trunc,100); 
  F2 = repmat(KL.mean,1,100) + KL.bases*diag(KL.sv)*W;
  %olot result
  figure(1)
  contourf(F); 
  
  % generate 2-D random field
  % build the correlation struct
  corr.name = 'exp';
  corr.c0 = [0.2 1]; % anisotropic correlation

  x = linspace(-1,1,11);
  [X,Y] = meshgrid(x,x); mesh = [X(:) Y(:)]; % 2-D mesh

  % set a spatially varying variance (must be positive!)
  corr.sigma = cos(pi*mesh(:,1)).*sin(2*pi*mesh(:,2))+1.5;

  [F,KL] = randomfield(corr,mesh,...
              'trunc', 10);

  % plot the realization
  figure(2)
  surf(X,Y,reshape(F,11,11)); view(2); colorbar;
  
  %random field fitting with kriging
     % create random field with autocorrelation
    [X,Y] = meshgrid(0:500);
    Z = randn(size(X));
    Z = imfilter(Z,fspecial('gaussian',[40 40],8));

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




