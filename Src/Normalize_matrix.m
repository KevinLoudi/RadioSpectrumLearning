% ****************************************************
% Note: normalize matrix data to [low, high]
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************
function mat_nor=Normalize_matrix(mat, low, high)
   %exception handling
   if nargin==0
       error('No input data!!');
   elseif nargin<3
      low=0.0; high=1.0;        
   end
   if (low>high)
       tmp=low; low=high; high=tmp;
   elseif (low==high)
       error('Two boundry cannot be assigned the same!!!');
   end
   
   %get matrix size
   [rows,cols]=size(mat);
   %calculate the max and min value
   mat_m=reshape(mat,rows*cols,1);
   max_v=max(mat_m);
   min_v=min(mat_m);
   %normailize
   mat_nor=((mat-min_v)/(max_v-min_v))*(high-low)+low;
end