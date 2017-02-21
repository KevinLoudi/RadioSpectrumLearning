%%
clear all;close all
TI=channels;
TI=TI(2:2:end,2:2:end);
SIM=ones(80,80)*NaN;

options.n_mulgrids=3;;
options.n_cond=[25 25 25 25];;
options.n_template=[9 9 9 9].^2;;
options.plot_interval=1e+9;;

rng(1);
[out_1,O_1]=mps_snesim(TI,SIM,options);


%%
O_2=options;
O_2.n_template=[9 9 9 9].^2;;
O_2.n_cond=[25 15 10 9];;
O_2.n_cond=[25 18 15 9];;

%O_2.T=O_1.T;
%O_2.ST_mul=O_1.ST_mul;
rng(1);[out_2,O_2]=mps_snesim(TI,SIM,O_2);

%%
figure(2);
subplot(1,3,1);
imagesc([out_1]);axis image;caxis([-1 1])
subplot(1,3,2);
imagesc([out_2]);axis image;;caxis([-1 1])
subplot(1,3,3);
imagesc([out_1-out_2]);axis image;;caxis([-1 1])

return
%%
rng(1);
clear all;close all
options.n_cond=55;
options.n_template=81;

%options.n_cond=19;
%options.n_template=19;

options.n_cond=49;
options.n_template=9;

%options.n_mulgrids=4;
options.n_mulgrids=4;

options.plot=1;options.plot_interval=100;

options.rand_path=1;

SIM=ones(200,400).*NaN;
TI=channels;TI=TI(4:4:end,4:4:end);
TI=channels;%TI=TI(2:2:end,2:2:end);

[out,o]=mps_snesim(TI,SIM,options);
return

figure(5);
subplot(2,3,1);imagesc(out);axis image
subplot(2,3,2);imagesc(o.C);colorbar;axis image
subplot(2,3,3);imagesc(o.IPATH);colorbar;axis image
options.n_max_ite=10000;
[out_dsim,o_dsim]=mps_enesim(TI,SIM,options);
subplot(2,3,4);imagesc(out_dsim);colorbar;axis image


return
%%
[out1,o1]=mps_snesim(TI,SIM,o);

%%
n=prod(size(SIM));
n_resim=ceil(0.05*n);
for i=1:50;
  SIM2=out1;
  
  i_resim=randomsample(1:n,n_resim);
  SIM2(i_resim)=NaN;
  [out1,o1]=mps_snesim(TI,SIM2,o1);
  figure(9)
  imagesc(out1);axis image;drawnow;
end
  



