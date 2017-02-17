clc
clear all

v = ver;
test_stat=any(strcmp('Statistics and Machine Learning Toolbox', {v.Name}));
if not(test_stat)
    error('You need the Statistics Toolbox and the Optimization Toolbox to run this case study. The Statistics Toolbox seems to be missing.');
end
test_opt=any(strcmp('Optimization Toolbox', {v.Name}));
if not(test_opt)
    error('You need the Statistics Toolbox and the Optimization Toolbox to run this case study. The Optimization Toolbox seems to be missing.');
end
test_map=any(strcmp('Mapping Toolbox', {v.Name}));
if not(test_map)
    disp('The Mapping Toolbox is missing. The output maps will not be displayed.');
    choice  =  input('Press any key to continue...');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Ground level data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the no2 ground level observations
load ../Data/no2_ground.mat
ground.Y{1} = no2_ground.data;
ground.Y_name{1} = 'no2 ground';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load the pm25 ground level observations
load ../Data/pm25_ground.mat
ground.Y{2} = pm25_ground.data;
ground.Y_name{2} = 'pm2.5 ground';
n2 = size(ground.Y{2},1);

%downscaler loading coefficients
ground.X_bp{1} = ones(n1, 1);
ground.X_bp_name{1} = {'constant'};
ground.X_bp{2} = ones(n2, 1);
ground.X_bp_name{2} = {'constant'};

%X_beta
%load the covariates for the NO2 monitoring stations
load ../Data/no2_ground_covariates.mat
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'wind speed', 'elevation', 'sunday'};

%load the covariates for the PM2.5 monitoring stations
load ../Data/pm25_ground_covariates.mat
ground.X_beta{2} = X;
ground.X_beta_name{2} = {'wind speed', 'elevation', 'sunday'};

%X_z
ground.X_z{1} = ones(n1, 1, 1);
ground.X_z_name{1} = {'constant'};
ground.X_z{2} = ones(n2, 1, 1);
ground.X_z_name{2} = {'constant'};

ground.X_p{1} = ground.X_beta{1}(:, 2, 1);
ground.X_p_name{1} = {'elevation'};
ground.X_p{2} = ground.X_beta{2}(:, 2, 1);
ground.X_p_name{2} = {'elevation'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, ground.X_bp, ...
                                ground.X_bp_name,ground.X_beta, ...
                                ground.X_beta_name,ground.X_z, ...
                                ground.X_z_name, ground.X_p, ...
                                ground.X_p_name);

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [no2_ground.lat, no2_ground.lon];
ground.coordinates{2} = [pm25_ground.lat, pm25_ground.lon];
obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Remote sensing data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Y
load ../Data/no2_remote_025.mat
remote.Y{1} = no2_remote.data;
remote.Y_name{1} = 'no2 remote';
m1 = size(remote.Y{1}, 1);

load ../Data/aot_remote_025.mat
remote.Y{2} = aot_remote.data;
remote.Y_name{2} = 'aot remote';
m2 = size(remote.Y{2}, 1);
    
%X_bp
remote.X_bp{1} = ones(m1, 1);
remote.X_bp_name{1} = {'constant'};
remote.X_bp{2} = ones(m2, 1);
remote.X_bp_name{2} = {'constant'};

obj_stem_varset_b = stem_varset(remote.Y, remote.Y_name, ...
                                remote.X_bp, remote.X_bp_name);

obj_stem_gridlist_b = stem_gridlist();
remote.coordinates{1} = [no2_remote.lat(:) ,no2_remote.lon(:)];
remote.coordinates{2} = [aot_remote.lat(:), aot_remote.lon(:)];
obj_stem_grid1 = stem_grid(remote.coordinates{1}, 'deg', 'regular', ...
                           'pixel', size(no2_remote.lat), 'square', ...
                           0.25, 0.25);
obj_stem_grid2 = stem_grid(remote.coordinates{2}, 'deg', 'regular', ...
                           'pixel', size(aot_remote.lat), 'square', ...
                           0.25, 0.25);
obj_stem_gridlist_b.add(obj_stem_grid1);
obj_stem_gridlist_b.add(obj_stem_grid2);
clear aot_remote
clear no2_remote

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2009 00:00', '31-12-2009 00:00', T);

flag_pixel_correlated = 0;
flag_time_diagonal = 0;

%stem_data object creation
if test_map
    shape = shaperead('../Maps/worldmap');
end
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          obj_stem_varset_b, obj_stem_gridlist_b, ...
                          obj_stem_datestamp, shape, [], [], ...
                          flag_pixel_correlated);
%stem_par object creation
obj_stem_par = stem_par(obj_stem_data, 'exponential', flag_time_diagonal);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
clear ground
clear remote

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.alpha_bp = [0.4 0.4 0.8 0.8]';
if flag_pixel_correlated
    obj_stem_par.theta_b = 100;
    obj_stem_par.v_b = [1 0.6; 0.6 1];
else
    obj_stem_par.theta_b = [100 100];
    obj_stem_par.v_b = eye(2);
end

obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.alpha_p = [0.6 0.6]';
obj_stem_par.theta_p = 100;
obj_stem_par.v_p = [1 0.6; 0.6 1];
obj_stem_par.sigma_eta = diag([0.2 0.2]);
obj_stem_par.G = diag([0.8 0.8]);
obj_stem_par.sigma_eps = diag([0.3 0.3 0.3 0.3]);

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 1;
obj_stem_EM_options = stem_EM_options(exit_toll, max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

obj_stem_model.print;        

if test_map
    obj_stem_model.stem_data.plot('no2 remote', 'pixel', 25);
    
    size = obj_stem_model.stem_data.stem_gridlist_b.grid{1}.grid_size;
    E_wb_y1 = obj_stem_model.stem_EM_result.E_wb_y1(1:size(1) * size(2), 25);
    coordinate = obj_stem_model.stem_data.stem_gridlist_b.grid{1}.coordinate;
    lat = coordinate(:, 1);
    lon = coordinate(:, 2);
    lat = reshape(lat, size);
    lon = reshape(lon, size);
    E_wb_y1 = reshape(E_wb_y1, size);
    stem_misc.plot_map(lat, lon, E_wb_y1, obj_stem_model.stem_data.shape, ...
                       'no2 remote estimated on 25-Jan-2009', ...
                       'Longitude', 'Latitude');
end

obj_stem_model_s6 = obj_stem_model;
save('../Output/obj_stem_model_s6.mat', 'obj_stem_model_s6', '-v7.3');