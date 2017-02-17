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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
n2 = size(ground.Y{2}, 1);

%load covariates for the NO2 monitoring stations
load ../Data/no2_ground_covariates.mat
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'wind speed', 'elevation', 'sunday'};

%load covariates for the PM2.5 monitoring stations
load ../Data/pm25_ground_covariates.mat
ground.X_beta{2} = X;
ground.X_beta_name{2} = {'wind speed', 'elevation', 'sunday'};

%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};
ground.X_z{2} = ones(n2, 1);
ground.X_z_name{2} = {'constant'};

ground.X_p{1} = ground.X_beta{1}(:, 2, 1);
ground.X_p_name{1} = {'elevation'};

ground.X_p{2} = ground.X_beta{2}(:, 2, 1);
ground.X_p_name{2} = {'elevation'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name, ...
                                ground.X_p, ground.X_p_name);

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [no2_ground.lat, no2_ground.lon];
ground.coordinates{2} = [pm25_ground.lat, pm25_ground.lon];
obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2009 00:00', '31-12-2009 00:00', T);

%stem_data object creation
if test_map
    shape = shaperead('../Maps/worldmap');
else
    shape=[];
end
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, shape);
%stem_par object creation
flag_time_diagonal = 0;
obj_stem_par = stem_par(obj_stem_data, 'exponential', flag_time_diagonal);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
clear ground

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.alpha_p = [0.6 0.6]';
obj_stem_par.theta_p = 100;
obj_stem_par.v_p = [1 0.6; 0.6 1];
obj_stem_par.sigma_eta = diag([0.2 0.2]);
obj_stem_par.G = diag([0.8 0.8]);
obj_stem_par.sigma_eps = diag([0.3 0.3]);

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 100;
obj_stem_EM_options = stem_EM_options(exit_toll, max_iterations, ...
                                      'single', [], 1);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;

obj_stem_model.print;        
obj_stem_model.stem_EM_result.stem_kalmansmoother_result.plot;

obj_stem_model_s5 = obj_stem_model;
save('../Output/obj_stem_model_s5.mat', 'obj_stem_model_s5', '-v7.3');