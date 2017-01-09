%no error handle at first
function mat_nor=normalize_matrix(mat, low, high)
   [rows,cols]=size(mat);
   mat_m=reshape(mat,rows*cols,1);
   max_v=max(mat_m);
   min_v=min(mat_m);
   mat_nor=((mat-min_v)/(max_v-min_v))*(high-low)+low;
end