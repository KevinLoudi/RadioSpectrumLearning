%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Economics and Quantitative Methods %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%% Code website: https://code.google.com/p/d-stem/                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.


classdef stem_EM < EM
    
    %CONSTANTS
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %S   = 2 if both point and pixel data are considered. S = 1 if only point data are considered.
    %T   - number of temporal steps
    %TT  = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    %p   - dimension of the latent temporal variable z
    
    properties
        stem_model=[];               %[stem_model object] (1x1)
        stem_EM_options=[];          %[stem_EM_options]   (1x1)
    end
    
    methods
        function obj = stem_EM(stem_model,stem_EM_options)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_model        - [stem_model object] (1x1)
            %stem_EM_options   - [stem_EM_options]   (1x1)
            %
            %OUTPUT
            %obj               - [stem_EM object]    (1x1)
            if nargin<2
                error('All the arguments must be provided');
            end
            
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('The first argument must be of class stem_model');
            end
            
            if isa(stem_EM_options,'stem_EM_options');
                obj.stem_EM_options=stem_EM_options;
            else
                error('The second argument must be of class stem_EM_options');
            end
            
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
        end
        
        function st_EM_result = estimate(obj)
            %DESCRIPTION: EM estimation
            %
            %INPUT
            %obj            - [stem_EM object]      (1x1)
            %
            %OUTPUT
            %st_EM_result   - [st_EM_result object] (1x1)
            
            t1_full=clock;
            if isempty(obj.stem_model)&&(nargin==0)
                error('You have to set the stem_model property first');
            end
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
            delta=9999;
            delta_logL=9999;
            last_logL=0;
            last_stem_par=obj.stem_model.stem_par;
            iteration=0;
            st_EM_result=stem_EM_result();
            st_EM_result.max_iterations=obj.stem_EM_options.max_iterations;
            st_EM_result.exit_toll=obj.stem_EM_options.exit_toll;
            st_EM_result.machine=computer;
            st_EM_result.date_start=datestr(now);
            while (delta>obj.stem_EM_options.exit_toll)&&(delta_logL>obj.stem_EM_options.exit_toll)&&(iteration<obj.stem_EM_options.max_iterations)
                ct1=clock;
                iteration=iteration+1;
                disp('************************');
                disp(['Iteration ',num2str(iteration),' started...']);
                disp('************************');
                
                clear E_wb_y1
                clear sum_Var_wb_y1
                clear diag_Var_wb_y1
                clear cov_wb_z_y1
                clear E_wg_y1
                clear sum_Var_wp_y1
                clear diag_Var_wp_y1
                clear cov_wp_z_y1
                clear M_cov_wb_wp_y1
                clear cov_wpk_wph_y1
                clear diag_Var_e_y1
                clear E_e_y1
                clear sigma_eps
                clear Xbeta
                
                [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result] = obj.E_step();
                model_changed = obj.M_step(E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,iteration);
                
                if (model_changed==0)
                    if not(isempty(st_kalmansmoother_result))
                        if not(st_kalmansmoother_result.logL==0)
                            logL=st_kalmansmoother_result.logL;
                            st_EM_result.logL_all(iteration)=logL;
                            delta_logL=abs(logL-last_logL)/abs(logL);
                            last_logL=logL;
                            disp('****************');
                            disp( ['logL: ',num2str(logL)]);
                            disp(['logL relative delta: ',num2str(delta_logL)]);
                        else
                            delta_logL=9999;
                        end
                    else
                        delta_logL=9999;
                    end
                else
                    delta_logL=9999;
                end
                if (model_changed==0)
                    delta=norm(obj.stem_model.stem_par.vec()-last_stem_par.vec())/norm(last_stem_par.vec());
                else
                    delta=9999;
                end
                last_stem_par=obj.stem_model.stem_par;
                disp(['Parameter delta norm: ',num2str(delta)]);
                obj.stem_model.stem_par.print;
                ct2=clock;
                disp('**********************************************');
                disp(['Iteration ',num2str(iteration),' ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                disp('**********************************************');
            end
            t2_full=clock;
            st_EM_result.stem_par=obj.stem_model.stem_par;
            st_EM_result.stem_kalmansmoother_result=st_kalmansmoother_result;
            st_EM_result.E_wp_y1=E_wp_y1;
            st_EM_result.E_wb_y1=E_wb_y1;
            st_EM_result.diag_Var_wp_y1=diag_Var_wp_y1;
            st_EM_result.diag_Var_wb_y1=diag_Var_wb_y1;
            st_EM_result.diag_Var_e_y1=diag_Var_e_y1;
            st_EM_result.y_hat=obj.stem_model.stem_data.Y;
            st_EM_result.y_hat(isnan(st_EM_result.y_hat))=0;
            st_EM_result.y_hat=st_EM_result.y_hat-E_e_y1;
            st_EM_result.res=obj.stem_model.stem_data.Y-st_EM_result.y_hat;
            st_EM_result.diag_Var_y_hat_back=zeros(size(st_EM_result.y_hat));

            blocks=[0 cumsum(obj.stem_model.dim)];
            counter=1;
            for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                if obj.stem_model.stem_data.stem_varset_p.standardized
                    s=obj.stem_model.stem_data.stem_varset_p.Y_stds{i};
                    m=obj.stem_model.stem_data.stem_varset_p.Y_means{i};
                end
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&not(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m;
                    y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                    var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                    var_y_hat_back=var_y_hat*s^2;
                end
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                    var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                    y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                    y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                    var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m)+(var_y_hat*s^2));
                end

                st_EM_result.y_hat_back(blocks(counter)+1:blocks(counter+1),:)=y_hat_back;
                st_EM_result.y_back(blocks(counter)+1:blocks(counter+1),:)=y_back;
                st_EM_result.res_back(blocks(counter)+1:blocks(counter+1),:)=y_back-y_hat_back;
                if not(isempty(var_y_hat_back))
                    st_EM_result.diag_Var_y_hat_back(blocks(counter)+1:blocks(counter+1),:)=var_y_hat_back;
                end
                counter=counter+1;
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    if obj.stem_model.stem_data.stem_varset_b.standardized
                        s=obj.stem_model.stem_data.stem_varset_b.Y_stds{i};
                        m=obj.stem_model.stem_data.stem_varset_b.Y_means{i};
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&not(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m;
                        y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                        var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                        var_y_hat_back=var_y_hat*s^2;
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                        y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                        y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                        y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                        var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m)+(var_y_hat*s^2));
                    end
                    st_EM_result.y_hat_back(blocks(counter)+1:blocks(counter+1),:)=y_hat_back;
                    st_EM_result.y_back(blocks(counter)+1:blocks(counter+1),:)=y_back;
                    if not(isempty(var_y_hat_back))
                        st_EM_result.diag_Var_y_hat_back(blocks(counter)+1:blocks(counter+1),:)=var_y_hat_back;
                    end
                    st_EM_result.res_back(blocks(counter)+1:blocks(counter+1),:)=y_back-y_hat_back;
                    counter=counter+1;
                end
            end
            
            st_EM_result.iterations=iteration;
            st_EM_result.computation_time=etime(t2_full,t1_full);
        end
        
        function st_EM_result = estimate_parallel(obj,path_distributed_computing)
            %DESCRIPTION: EM parallel estimation
            %
            %INPUT
            %obj                                - [stem_EM object]      (1x1)
            %path_distributed_computing         - [string]              (1x1)  full or relative path of the folder to use for distributed computing
            %
            %OUTPUT
            %st_EM_result         - [st_EM_result object] (1x1)
            t1_full=clock;
            if isempty(obj.stem_model)&&(nargin==0)
                error('You have to set the stem_model property first');
            end
            if isempty(obj.stem_model.stem_par_initial)
                error('Initial value estimation for model parameters must be provided first');
            end
            
            T=obj.stem_model.T;
            K=obj.stem_model.stem_par.k;
            local_efficiency=1;
            delta=9999;
            delta_logL=9999;
            last_logL=0;
            last_stem_par=obj.stem_model.stem_par;
            iteration=0;
            st_EM_result=stem_EM_result();
            st_EM_result.max_iterations=obj.stem_EM_options.max_iterations;
            st_EM_result.exit_toll=obj.stem_EM_options.exit_toll;
            st_EM_result.machine=computer;
            st_EM_result.date_start=datestr(now);
            while (delta>obj.stem_EM_options.exit_toll)&&(delta_logL>obj.stem_EM_options.exit_toll)&&(iteration<obj.stem_EM_options.max_iterations)
                ct1_iteration=clock;
                iteration=iteration+1;
                disp('************************');
                disp(['Iteration ',num2str(iteration),' started...']);
                disp('************************');
                
                %repeat the E-step until no timeout occurs
                timeout=1;
                while timeout
                    %set the timeout to zero
                    timeout=0;
                    clear E_wb_y1
                    clear sum_Var_wb_y1
                    clear diag_Var_wb_y1
                    clear cov_wb_z_y1
                    clear E_wp_y1
                    clear sum_Var_wp_y1
                    clear diag_Var_wp_y1
                    clear cov_wp_z_y1
                    clear M_cov_wb_wp_y1
                    clear cov_wpk_wph_y1
                    clear diag_Var_e_y1
                    clear E_e_y1
                    clear sigma_eps
                    
                    %delete all the file in the exchange directory
                    files=dir([path_distributed_computing,'*.mat']);
                    for i=1:length(files)
                        delete([path_distributed_computing,files(i).name]);
                    end
                    
                    %create the file for the whoishere request
                    disp('  Looking for slaves...');
                    whoishere.IDrequest=unifrnd(0,100000,1,1);
                    
                    save([path_distributed_computing,'temp/whoishere.mat'],'whoishere');
                    pause(0.5);
                    movefile([path_distributed_computing,'temp/whoishere.mat'],[path_distributed_computing,'whoishere.mat']);
                    
                    if iteration==1
                        %hosts=[];
                        hosts=struct([]);
                    end
                    nhosts=length(hosts);
                    
                    %wait for the replies from the slaves
                    wait1=clock;
                    exit=0;
                    while not(exit)
                        files=dir([path_distributed_computing,'machine_*.*']);
                        for i=1:length(files)
                            try
                                load([path_distributed_computing,files(i).name])
                                if machine.IDrequest==whoishere.IDrequest
                                    %check if the slave is already in the hosts list
                                    idx=[];
                                    for j=1:nhosts
                                        if hosts(j).node_code==machine.node_code
                                            hosts(j).active=1;
                                            hosts(j).require_stemmodel=machine.require_stemmodel;
                                            idx=j;
                                        end
                                    end
                                    %if not, add the slave
                                    if isempty(idx)
                                        nhosts=nhosts+1;
                                        hosts(nhosts).node_code=machine.node_code;
                                        hosts(nhosts).data_received=0;
                                        %the first time a slave is added it has the efficiency of the server
                                        hosts(nhosts).efficiency=local_efficiency;
                                        hosts(nhosts).require_stemmodel=machine.require_stemmodel;
                                        hosts(nhosts).active=1;
                                    end
                                end
                            catch
                            end
                        end
                        wait2=clock;
                        if etime(wait2,wait1)>obj.stem_EM_options.timeout_node_search
                            exit=1;
                        end
                        pause(0.1);
                    end
                    delete([path_distributed_computing,'whoishere.mat']);
                    %check for inactive slaves
                    idx=[];
                    for i=1:nhosts
                        if hosts(i).active==0
                            idx=cat(2,idx,i);
                        end
                    end
                    if not(isempty(idx))
                        hosts(idx)=[];
                        nhosts=length(hosts);
                    end
                    
                    %if there is at least one slave then distribute the st_model
                    if nhosts>=1
                        disp(['  ',num2str(nhosts),' slave(s) found']);
                        disp('  Sending stem_model object to slave(s)...');
                        st_model=obj.stem_model;
                        compute_logL=obj.stem_EM_options.compute_logL_at_all_steps;
                        for i=1:nhosts
                            if hosts(i).require_stemmodel
                                save([path_distributed_computing,'temp/st_model_parallel_',num2str(hosts(i).node_code),'.mat'],'st_model','compute_logL','-v7.3');
                                pause(0.5);
                                movefile([path_distributed_computing,'temp/st_model_parallel_',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'st_model_parallel_',num2str(hosts(i).node_code),'.mat']);
                            end
                        end
                    else
                        disp('  No slaves found. Local computing.');
                    end
                    
                    if nhosts>=1
                        %the st_par to be distributed is the same for all
                        %the slaves
                        disp('  Sending stem_par object to slave(s)...')
                        st_par=obj.stem_model.stem_par;
                        for i=1:nhosts
                            save([path_distributed_computing,'temp/st_par_parallel_',num2str(hosts(i).node_code),'.mat'],'st_par');
                            pause(0.5);
                            movefile([path_distributed_computing,'temp/st_par_parallel_',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'st_par_parallel_',num2str(hosts(i).node_code),'.mat']);
                        end
                        clear st_par
                    end
                    
                    Lt_all=sum(not(isnan(obj.stem_model.stem_data.Y)));
                    Lt_sum=sum(Lt_all);
                    Lt_csum=cumsum(Lt_all);
                    veff=local_efficiency;
                    for i=1:nhosts
                        veff=cat(2,veff,hosts(i).efficiency);
                    end
                    veff=veff/sum(veff);
                    veff=[0 cumsum(veff)];
                    if not(veff(end)==1)
                        veff(end)=1;
                    end
                    %compute the time_steps for the server
                    l1=Lt_sum*veff(1);
                    l2=Lt_sum*veff(2);
                    t1=find(Lt_csum>l1,1);
                    t2=find(Lt_csum>=l2,1);
                    time_steps=t1:t2;
                    local_cb=sum(Lt_all(time_steps));
                    disp(['  ',num2str(length(time_steps)),' time steps will be assigned to the master node']);
                    
                    %Kalman smoother
                    if obj.stem_model.stem_par.p>0
                        %distribute the st_par and the data needed to the
                        %slaves
                        if nhosts>=1
                            disp('  Sending Kalman filter data to slave(s)...')
                            %send the information for the computation of the parallel kalman
                            data.iteration=iteration;
                            last_t2=t2;
                            for i=1:nhosts
                                %compute the time_steps for the slaves
                                l1=Lt_sum*veff(i+1);
                                l2=Lt_sum*veff(i+2);
                                t1=find(Lt_csum>l1,1);
                                if t1<=last_t2
                                    t1=last_t2+1;
                                end
                                t2=find(Lt_csum>=l2,1);
                                if t2<t1
                                    t2=t1;
                                end
                                data.time_steps=t1:t2;
                                last_t2=t2;
                                disp(['  ',num2str(length(data.time_steps)),' time steps assigned to slave ',num2str(hosts(i).node_code)]);
                                save([path_distributed_computing,'temp/kalman_parallel_',num2str(hosts(i).node_code),'.mat'],'data');
                                pause(0.5);
                                movefile([path_distributed_computing,'temp/kalman_parallel_',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'kalman_parallel_',num2str(hosts(i).node_code),'.mat']);
                            end
                            %local Kalman Smoother computation
                            st_kalman=stem_kalman(obj.stem_model);
                            [st_kalmansmoother_result,sigma_eps,~,~,~,~,~,~,~] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps,0,time_steps,path_distributed_computing);
                        else
                            %The computation is only local. The standard Kalman smoother is considered
                            st_kalman=stem_kalman(obj.stem_model);
                            [st_kalmansmoother_result,sigma_eps,~,~,~,~,~,~,~] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps,0);
                            time_steps=1:T;
                        end
                    else
                        st_kalmansmoother_result=[];
                        %sigma_eps
                        d=[];
                        dim=obj.stem_model.stem_data.dim;
                        for i=1:obj.stem_model.stem_data.nvar
                            d=cat(1,d,repmat(obj.stem_model.stem_par.sigma_eps(i,i),dim(i),1));
                        end
                        sigma_eps=diag(d);
                    end
                    
                    disp('  Sending E-step data file to slave(s)...')
                    data.st_kalmansmoother_result=st_kalmansmoother_result;
                    data.iteration=iteration;
                    
                    l2=Lt_sum*veff(2);
                    t2=find(Lt_csum>=l2,1);
                    last_t2=t2;
                    for i=1:nhosts
                        %compute the time_steps for the slaves
                        l1=Lt_sum*veff(i+1);
                        l2=Lt_sum*veff(i+2);
                        t1=find(Lt_csum>l1,1);
                        if t1<=last_t2
                            t1=last_t2+1;
                        end
                        t2=find(Lt_csum>=l2,1);
                        if t2<t1
                            t2=t1;
                        end
                        data.time_steps=t1:t2;
                        data.cb=sum(Lt_all(data.time_steps));
                        last_t2=t2;
                        disp(['  ',num2str(length(data.time_steps)),' time steps assigned to slave ',num2str(hosts(i).node_code)]);
                        save([path_distributed_computing,'temp/data_parallel_',num2str(hosts(i).node_code),'.mat'],'data');
                        pause(0.5);
                        movefile([path_distributed_computing,'temp/data_parallel_',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'data_parallel_',num2str(hosts(i).node_code),'.mat']);
                    end
                    clear data
                    
                    %local E-step computation
                    ct1_local=clock;
                    [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1] = obj.E_step_parallel(time_steps,st_kalmansmoother_result);
                    ct2_local=clock;
                    local_efficiency=local_cb/etime(ct2_local,ct1_local);
                    disp(['    Local computation time: ',num2str(etime(ct2_local,ct1_local))]);
                    disp(['    Local efficiency: ',num2str(local_efficiency)]);
                    
                    if nhosts>=1
                        disp('    Waiting for the results from slave(s)...');
                        exit=0;
                        wait1=clock;
                        while not(exit)
                            files=dir([path_distributed_computing,'output_*.*']);
                            for i=1:length(files)
                                load([path_distributed_computing,files(i).name]);
                                disp(['    Received result from slave ',num2str(output.node_code)]);
                                if iteration==output.iteration
                                    idx=[];
                                    for j=1:nhosts
                                        if (hosts(j).node_code==output.node_code)&&(hosts(j).data_received==0)
                                            idx=j;
                                        end
                                    end
                                    if not(isempty(idx))
                                        disp('    The result data file from the slave was expected.');
                                        hosts(idx).efficiency=output.cb/output.ct;
                                        disp(['    Computational time of slave ',num2str(hosts(idx).node_code),': ',num2str(output.ct)]);
                                        disp(['    Efficiency of slave ',num2str(hosts(idx).node_code),': ',num2str(hosts(idx).efficiency)]);
                                        tsteps=output.time_steps;
                                        if not(isempty(E_wb_y1))
                                            E_wb_y1(:,tsteps)=output.E_wb_y1;
                                        end
                                        if not(isempty(sum_Var_wb_y1))
                                            %the matrix is recomposed since only the upper triangular part is received
                                            sum_Var_wb_y1=sum_Var_wb_y1+output.sum_Var_wb_y1+triu(output.sum_Var_wb_y1,1)';
                                        end
                                        if not(isempty(diag_Var_wb_y1))
                                            diag_Var_wb_y1(:,tsteps)=output.diag_Var_wb_y1;
                                        end
                                        if not(isempty(cov_wb_z_y1))
                                            cov_wb_z_y1(:,:,tsteps)=output.cov_wb_z_y1;
                                        end
                                        if not(isempty(E_wp_y1))
                                            for k=1:K
                                                E_wp_y1(:,tsteps,k)=output.E_wp_y1(:,:,k);
                                            end
                                        end
                                        if not(isempty(sum_Var_wp_y1))
                                            for k=1:K
                                                %the matrix is recomposed since only the upper triangular part is received
                                                sum_Var_wp_y1{k}=sum_Var_wp_y1{k}+output.sum_Var_wp_y1{k}+triu(output.sum_Var_wp_y1{k},1)';
                                            end
                                        end
                                        if not(isempty(diag_Var_wp_y1))
                                            for k=1:K
                                                diag_Var_wp_y1(:,tsteps,k)=output.diag_Var_wp_y1(:,:,k);
                                            end
                                        end
                                        if not(isempty(cov_wp_z_y1))
                                            for k=1:K
                                                cov_wp_z_y1(:,:,tsteps,k)=output.cov_wp_z_y1(:,:,:,k);
                                            end
                                        end
                                        if not(isempty(M_cov_wb_wp_y1))
                                            for k=1:K
                                                M_cov_wb_wp_y1(:,tsteps,k)=output.M_cov_wb_wp_y1(:,:,k);
                                            end
                                        end
                                        if iscell(cov_wpk_wph_y1)
                                            for h=1:K
                                                for k=h+1:K
                                                    cov_wpk_wph_y1{k,h}(:,tsteps)=output.cov_wpk_wph_y1{k,h};
                                                end
                                            end
                                        end
                                        diag_Var_e_y1(:,tsteps)=output.diag_Var_e_y1;
                                        E_e_y1(:,tsteps)=output.E_e_y1;
                                        
                                        hosts(idx).data_received=1;
                                        clear output
                                    else
                                        disp('    Something is wbong');
                                    end
                                    exit=1;
                                    for j=1:nhosts
                                        if hosts(j).data_received==0
                                            exit=0;
                                        end
                                    end
                                    if exit==1
                                        disp('    All the data from the slave(s) have been collected');
                                    end
                                else
                                    disp('    The iteration within the output file does not match. The file is deleted');
                                end
                                deleted=0;
                                while not(deleted)
                                    try
                                        delete([path_distributed_computing,files(i).name]);
                                        deleted=1;
                                    catch
                                    end
                                end
                            end
                            wait2=clock;
                            if etime(wait2,wait1)>obj.stem_EM_options.timeout_distributed_computing
                                disp('    Timeout');
                                timeout=1;
                                exit=1;
                            end
                            pause(0.02);
                        end
                    end
                    
                    for i=1:nhosts
                        hosts(i).active=0;
                        hosts(i).data_received=0;
                    end
                end
                
                clear data
                if (K<=1)
                    %run the non parallel version of the M-step
                    model_changed = obj.M_step(E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,iteration);
                    %send the message to the other machine that they don't have to run the M-step
                    for i=1:nhosts
                        data.iteration=iteration;
                        data.index=[];
                        save([path_distributed_computing,'temp/data_parallel_mstep',num2str(hosts(i).node_code),'.mat'],'data');
                        pause(0.5);
                        movefile([path_distributed_computing,'temp/data_parallel_mstep',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'data_parallel_mstep',num2str(hosts(i).node_code),'.mat']);
                    end
                else
                    step=ceil(K/(nhosts+1));
                    counter=1;
                    index=cell(nhostrs,1);
                    for i=1:nhosts+1
                        if (i==1)
                            index_local=counter:step+counter-1;
                        else
                            index{i-1}=counter:step+counter-1;
                            index{i-1}(index{i-1}>K)=[];
                        end
                        counter=counter+step;
                    end
                    %send the messages to the host
                    clear data
                    for i=1:nhosts
                        data.iteration=iteration;
                        data.index=index{i};
                        disp(['    Preparing M-step data file for slave ',num2str(hosts(i).node_code)]);
                        data.sum_Var_wp_y1=sum_Var_wp_y1(index{i});
                        data.E_wp_y1=E_wp_y1(:,:,index{i});
                        disp(['    Sending M-step data to slave ',num2str(hosts(i).node_code)]);
                        save([path_distributed_computing,'temp/data_parallel_mstep',num2str(hosts(i).node_code),'.mat'],'data');
                        pause(0.5);
                        movefile([path_distributed_computing,'temp/data_parallel_mstep',num2str(hosts(i).node_code),'.mat'],[path_distributed_computing,'data_parallel_mstep',num2str(hosts(i).node_code),'.mat']);
                        disp('    M-Step data sent.');
                    end
                    %M-step locale
                    model_changed = obj.M_step_parallel(E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,index_local);
                end
                
                %Wait for the other nodes
                if nhosts>0
                    disp('  Waiting for M-step result files from the slave(s)...');
                    exit=0;
                    while not(exit)
                        files=dir([path_distributed_computing,'output_mstep_*.*']);
                        for i=1:length(files)
                            % try
                            load([path_distributed_computing,files(i).name]);
                            disp(['  Received M-step result file from slave ',num2str(output.node_code)]);
                            if iteration==output.iteration
                                idx=[];
                                for j=1:nhosts
                                    if (hosts(j).node_code==output.node_code)&&(hosts(j).data_received==0)
                                        idx=j;
                                    end
                                end
                                if not(isempty(idx))
                                    disp('  The M-step result file from the slave was expected');
                                    if not(isempty(output.index))
                                        for z=1:length(output.index)
                                            obj.stem_model.stem_par.v_p(:,:,output.index(z))=output.mstep_par.v_p(:,:,output.index(z));
                                            obj.stem_model.stem_par.theta_p(:,output.index(z))=output.mstep_par.theta_p(:,output.index(z));
                                            disp(['  ',num2str(output.index(z)),'th component of vg and theta_p updated']);
                                        end
                                    else
                                        disp('  The M-step result file from the slave is empty as expected');
                                    end
                                    hosts(idx).data_received=1;
                                    clear output
                                else
                                    disp('    Something went wrong!');
                                end
                                exit=1;
                                for j=1:nhosts
                                    if hosts(j).data_received==0
                                        exit=0;
                                    end
                                end
                                if exit==1
                                    disp('  All the M-step result data files from the slave(s) have been collected');
                                end
                            else
                                disp('    The iteration number within the output file does not match. The file is deleted');
                            end
                            deleted=0;
                            while not(deleted)
                                try
                                    delete([path_distributed_computing,files(i).name]);
                                    deleted=1;
                                catch
                                end
                            end
                            %catch
                            %end
                        end
                        wait2=clock;
                        if etime(wait2,wait1)>obj.stem_EM_options.timeout_distributed_computing
                            disp('    Timeout');
                            exit=1;
                        end
                        pause(0.05);
                    end
                    for i=1:nhosts
                        hosts(i).active=0;
                        hosts(i).data_received=0;
                    end
                end
                
                if model_changed==0
                    if not(isempty(st_kalmansmoother_result))
                        if not(st_kalmansmoother_result.logL==0)
                            logL=st_kalmansmoother_result.logL;
                            st_EM_result.logL_all(iteration)=logL;
                            delta_logL=abs(logL-last_logL)/abs(logL);
                            last_logL=logL;
                            disp('****************');
                            disp( ['logL: ',num2str(logL)]);
                            disp(['logL relative delta: ',num2str(delta_logL)]);
                        else
                            delta_logL=9999;
                        end
                    else
                        delta_logL=9999;
                    end
                else
                    delta_logL=9999;
                end
                if model_changed==0
                    delta=norm(obj.stem_model.stem_par.vec()-last_stem_par.vec())/norm(last_stem_par.vec());
                else
                    delta=9999;
                end
                last_stem_par=obj.stem_model.stem_par;
                disp(['Parameter delta norm: ',num2str(delta)]);
                obj.stem_model.stem_par.print;
                ct2_iteration=clock;
                disp('**********************************************');
                disp(['Iteration ',num2str(iteration),' ended in ',stem_misc.decode_time(etime(ct2_iteration,ct1_iteration))]);
                disp('**********************************************');
            end
            t2_full=clock;
            
            st_EM_result.stem_par=obj.stem_model.stem_par;
            st_EM_result.stem_kalmansmoother_result=st_kalmansmoother_result;
            st_EM_result.E_wp_y1=E_wp_y1;
            st_EM_result.E_wb_y1=E_wb_y1;
            st_EM_result.diag_Var_wp_y1=diag_Var_wp_y1;
            st_EM_result.diag_Var_wb_y1=diag_Var_wb_y1;
            st_EM_result.diag_Var_e_y1=diag_Var_e_y1;
            st_EM_result.y_hat=obj.stem_model.stem_data.Y;
            st_EM_result.y_hat(isnan(st_EM_result.y_hat))=0;
            st_EM_result.y_hat=st_EM_result.y_hat-E_e_y1;
            st_EM_result.res=obj.stem_model.stem_data.Y-st_EM_result.y_hat;

            blocks=[0 cumsum(obj.stem_model.dim)];
            counter=1;
            for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                if obj.stem_model.stem_data.stem_varset_p.standardized
                    s=obj.stem_model.stem_data.stem_varset_p.Y_stds{i};
                    m=obj.stem_model.stem_data.stem_varset_p.Y_means{i};
                end
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&not(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m;
                    y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                    var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                    var_y_hat_back=var_y_hat*s^2;
                end
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                    y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                    y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                    y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                    var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m)+(var_y_hat*s^2));
                end

                st_EM_result.y_hat_back(blocks(counter)+1:blocks(counter+1),:)=y_hat_back;
                st_EM_result.y_back(blocks(counter)+1:blocks(counter+1),:)=y_back;
                if not(isempty(var_y_hat_back))
                    st_EM_result.diag_Var_y_hat_back(blocks(counter)+1:blocks(counter+1),:)=var_y_hat_back;
                end
                st_EM_result.res_back(blocks(counter)+1:blocks(counter+1),:)=y_back-y_hat_back;
                counter=counter+1;
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:);
                var_y_hat_back=[];
                
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    if obj.stem_model.stem_data.stem_varset_b.standardized
                        s=obj.stem_model.stem_data.stem_varset_b.Y_stds{i};
                        m=obj.stem_model.stem_data.stem_varset_b.Y_means{i};
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&not(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m;
                        y_back=obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m;
                        var_y_hat=st_EM_result.diag_Var_e_y1(blocks(counter)+1:blocks(counter+1),:);
                        var_y_hat_back=var_y_hat*s^2;
                    end
                    if (obj.stem_model.stem_data.stem_varset_b.standardized)&&(obj.stem_model.stem_data.stem_varset_b.log_transformed)
                        y_hat_back=st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:);
                        y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                        y_back=exp(obj.stem_model.stem_data.Y(blocks(counter)+1:blocks(counter+1),:)*s+m);
                        var_y_hat_back=(exp(var_y_hat*s^2)-1).*exp(2*(st_EM_result.y_hat(blocks(counter)+1:blocks(counter+1),:)*s+m)+(var_y_hat*s^2));
                    end
                    st_EM_result.y_hat_back(blocks(counter)+1:blocks(counter+1),:)=y_hat_back;
                    st_EM_result.y_back(blocks(counter)+1:blocks(counter+1),:)=y_back;
                    if not(isempty(var_y_hat_back))
                        st_EM_result.diag_Var_y_hat_back(blocks(counter)+1:blocks(counter+1),:)=var_y_hat_back;
                    end
                    st_EM_result.res_back(blocks(counter)+1:blocks(counter+1),:)=y_back-y_hat_back;
                    counter=counter+1;
                end
            end
         
            st_EM_result.iterations=iteration;
            st_EM_result.computation_time=etime(t2_full,t1_full);
        end
        
        function [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result] = E_step(obj,T)
            %DESCRIPTION: E-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1)
            %<T>                            - [integer >0]      (1x1) The E-step is computed only for the data related to the time steps between 1 and T
            %
            %OUTPUT
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp_k|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            %sigma_eps                      - [double]          (NxN) sigma_eps
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)
            
            if nargin==1
                T=obj.stem_model.stem_data.T;
            end
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            par=obj.stem_model.stem_par;
            
            disp('  E step started...');
            ct1_estep=clock;
          
            if p>0
                %Kalman smoother
                st_kalman=stem_kalman(obj.stem_model);
                [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,aj_p,aj_z,M] = st_kalman.smoother(obj.stem_EM_options.compute_logL_at_all_steps,0);
                
                rr=size(sigma_Z,1);
                
                if not(obj.stem_model.stem_data.X_z_tv)
                    if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                        temp=obj.stem_model.stem_data.X_z(:,:,1);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,1);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,1),1),size(obj.stem_model.stem_data.X_z(:,:,1),2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    
                    if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                        if obj.stem_model.tapering
                            %is it possible to improve the sparse matrix var_Zt?
                            var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                        else
                            var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                        end
                    else
                        var_Zt=[];
                    end
                end
                if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                    var_Yt=sigma_geo+var_Zt;
                end
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,~,~,aj_bp,aj_p,~,M] = obj.stem_model.get_sigma();

                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[],[],[]);
                var_Zt=[];
                rr=0;
                
                if not(isempty(sigma_geo))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
            end
            
            if (obj.stem_model.stem_par.model_type==2||obj.stem_model.stem_par.model_type==3)
                obj.stem_model.stem_data.stem_varset_p.Y{1}=[];
            end 
            
            E_e_y1=obj.stem_model.stem_data.Y;
            E_e_y1(isnan(E_e_y1))=0;
            if not(isempty(obj.stem_model.stem_data.X_beta))
                disp('    Xbeta evaluation started...');
                ct1=clock;
                Xbeta=zeros(N,T);
                if obj.stem_model.stem_data.X_beta_tv
                    for t=1:T
                        if size(obj.stem_model.stem_data.X_beta(:,:,t),1)<N
                            X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,t);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,t),1),size(obj.stem_model.stem_data.X_beta(:,:,t),2))];
                        else
                            X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,t);
                        end
                        Xbeta(:,t)=X_beta_orlated*par.beta;
                    end
                else
                    if size(obj.stem_model.stem_data.X_beta(:,:,1),1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,1);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,1),1),size(obj.stem_model.stem_data.X_beta(:,:,1),2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,1);
                    end
                    Xbeta=repmat(X_beta_orlated*par.beta,1,T);
                end
                ct2=clock;
                disp(['    Xbeta evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                E_e_y1=E_e_y1-Xbeta;
            else
                Xbeta=[];
            end
            diag_Var_e_y1=zeros(N,T,'single');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %sigma_Z=Var(Zt)
            %var_Zt=Var(X_z*Zt*X_z')
            
            disp('    Conditional E, Var, Cov evaluation started...');
            ct1=clock;
            if not(isempty(obj.stem_model.stem_data.X_bp))
                if (obj.stem_model.tapering)
                    Lr=find(sigma_W_b);
                    nnz_b=length(Lr);
                end
                %cov_wb_yz time invariant case
                if not(obj.stem_model.stem_data.X_bp_tv)
                    cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp(:,1,1),'r'),aj_bp,'r');
                end
                E_wb_y1=zeros(Nb,T);
                if (obj.stem_model.tapering)
                    sum_Var_wb_y1=spalloc(size(sigma_W_b,1),size(sigma_W_b,2),nnz_b);
                else
                    sum_Var_wb_y1=zeros(Nb);
                end
                diag_Var_wb_y1=zeros(Nb,T);
                cov_wb_z_y1=zeros(Nb,rr,T);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_p))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_p{1});
                    nnz_p=length(Lg);
                end
                %cov_wp_yz time invariant case
                if not(obj.stem_model.stem_data.X_p_tv)
                    cov_wp_y=cell(K,1);
                    for k=1:K
                        cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,1,k),'r'),aj_p(:,k),'r');
                    end
                end
                cov_wpk_wph_y1=cell(K,K);
                for h=1:K
                    for k=h+1:K
                        %these matrices can be sparse?
                        cov_wpk_wph_y1{k,h}=zeros(Np,T);
                    end
                end
                E_wp_y1=zeros(Np,T,K);
                sum_Var_wp_y1=cell(K,1);
                for k=1:K
                    if obj.stem_model.tapering
                        sum_Var_wp_y1{k}=spalloc(size(sigma_W_p{k},1),size(sigma_W_p{k},2),nnz_p);
                    else
                        sum_Var_wp_y1{k}=zeros(Np,Np);
                    end
                end
                diag_Var_wp_y1=zeros(Np,T,K);
                cov_wp_z_y1=zeros(Np,rr,T,K);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_bp)) && not(isempty(obj.stem_model.stem_data.X_p))
                M_cov_wb_wp_y1=zeros(N,T,K);
            else
                M_cov_wb_wp_y1=[];
            end
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                
                if obj.stem_model.stem_data.X_bp_tv
                    tBP=t;
                else
                    tBP=1;
                end
                if obj.stem_model.stem_data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if obj.stem_model.stem_data.X_p_tv
                    tP=t;
                else
                    tP=1;
                end
                
                %evaluate var_yt in the time variant case
                if obj.stem_model.stem_data.X_tv
                    if not(isempty(obj.stem_model.stem_data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                    end
                    if not(isempty(obj.stem_model.stem_data.X_p))
                        if isempty(sigma_geo)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(obj.stem_model.stem_data.X_p,4)
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,tP,k),'b'),aj_p(:,k),'b');
                        end
                    end
                    if isempty(sigma_geo)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if p>0
                        if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                            temp=obj.stem_model.stem_data.X_z(:,:,tT);
                            temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                        end
                        X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                        
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            if obj.stem_model.tapering
                                %is it possible to improve the sparse matrix var_Zt?
                                var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                            else
                                var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                            end
                            if isempty(sigma_geo)
                                var_Yt=var_Zt;
                            else
                                var_Yt=sigma_geo+var_Zt;
                            end
                        else
                            var_Zt=[];
                            var_Yt=[];
                        end
                    else
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            var_Yt=sigma_geo;
                        else
                            var_Yt=[];
                        end
                    end
                end
                
                if p>0
                    if N>obj.stem_model.system_size
                        blocks=0:80:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(X_z_orlated(blocks(i)+1:blocks(i+1),:)*st_kalmansmoother_result.Pk_s(:,:,t+1)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                        end
                    else
                        temp=X_z_orlated*st_kalmansmoother_result.Pk_s(:,:,t+1);
                        diag_Var_e_y1(:,t)=diag(temp*X_z_orlated');
                    end
                    %update E(e|y1)
                    temp=st_kalmansmoother_result.zk_s(:,t+1);
                    E_e_y1(:,t)=E_e_y1(:,t)-X_z_orlated*temp;
                end
                
                if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), X_z_orlated(Lt,:)*sigma_Z; sigma_Z*X_z_orlated(Lt,:)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    
                    res=obj.stem_model.stem_data.Y;
                    if not(isempty(Xbeta))
                        res=res-Xbeta;
                    end
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=stem_misc.chol_solve(chol_H1t,temp2(r));
                    else
                        chol_H1t=chol(H1t);
                        cs=stem_misc.chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end

                if not(isempty(obj.stem_model.stem_data.X_bp))
                    %cov_wb_yz time variant case
                    if obj.stem_model.stem_data.X_bp_tv
                        cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp(:,1,tBP),'r'),aj_bp,'r');
                    end
                    cov_wb_y1z=[cov_wb_y(:,Lt),zeros(size(cov_wb_y,1),rr)];
                    %compute E(w_b|y1);
                    E_wb_y1(:,t)=cov_wb_y1z*cs;
                    %compute Var(w_b|y1)
                    if obj.stem_model.tapering
                        temp_b(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wb_y1z(:,r)');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    else
                        temp_b=stem_misc.chol_solve(chol_H1t,cov_wb_y1z');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    end
                    
                    if p>0
                        %compute cov(w_b,z|y1)
                        cov_wb_z_y1(:,:,t)=temp_b(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                        Var_wb_y1=Var_wb_y1+cov_wb_z_y1(:,:,t)*temp_b(end-rr+1:end,:);
                        %update diag(Var(e|y1))
                        temp=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                        if N>obj.stem_model.system_size
                            blocks=0:80:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                        end
                    else
                        cov_wb_z_y1=[];
                    end
                    %compute diag(Var(w_b|y1))
                    diag_Var_wb_y1(:,t)=diag(Var_wb_y1);
                    %compute sum(Var(w_b|y1))
                    sum_Var_wb_y1=sum_Var_wb_y1+Var_wb_y1;
                    %update E(e|y1)
                    E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(diag_Var_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                else
                    E_wb_y1=[];
                    diag_Var_wb_y1=[];
                    sum_Var_wb_y1=[];
                    cov_wb_z_y1=[];
                end
                clear temp_b
                
                if not(isempty(obj.stem_model.stem_data.X_p))
                    if obj.stem_model.stem_data.X_p_tv
                        %cov_wp_yz time invariant case
                        cov_wp_y=cell(K,1);
                        for k=1:K
                            cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,tP,k),'r'),aj_p(:,k),'r');
                        end
                    end
                    temp_p=cell(K,1);
                    for k=1:K
                        cov_wp_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                        %compute E(w_p_k|y1);
                        E_wp_y1(:,t,k)=cov_wp_y1z*cs;
                        %compute Var(w_p_k|y1)
                        if obj.stem_model.tapering
                            temp_p{k}(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wp_y1z(:,r)');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        else
                            temp_p{k}=stem_misc.chol_solve(chol_H1t,cov_wp_y1z');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        end
                        
                        if p>0
                            %compute cov(w_p,z|y1)
                            cov_wp_z_y1(:,:,t,k)=temp_p{k}(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                            Var_wp_y1=Var_wp_y1+cov_wp_z_y1(:,:,t,k)*temp_p{k}(end-rr+1:end,:);
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                            end
                        else
                            cov_wp_z_y1=[];
                        end
                        diag_Var_wp_y1(:,t,k)=diag(Var_wp_y1);
                        sum_Var_wp_y1{k}=sum_Var_wp_y1{k}+Var_wp_y1;
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(diag_Var_wp_y1(:,t,k),obj.stem_model.stem_data.X_p(:,:,tP,k),'b'),aj_p(:,k),'b'); %K varianze
                        
                        if not(isempty(obj.stem_model.stem_data.X_bp))
                            %compute M_cov(w_b,w_p|y1) namely M*cov(w_b,w_p|y1)
                            if length(M)>obj.stem_model.system_size
                                blocks=0:80:length(M);
                                if not(blocks(end)==length(M))
                                    blocks=cat(2,blocks,length(M));
                                end
                                for i=1:length(blocks)-1
                                    if p>0
                                        M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1))+cov_wb_z_y1(M(blocks(i)+1:blocks(i+1)),:,t)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1))); %ha gia' l'stem_misc.M_apply su left!!
                                    else
                                        M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M))+cov_wb_z_y1(M,:,t)*temp_p{k}(end-rr+1:end,1:length(M))); %ha gi� l'stem_misc.M_apply su left!!
                                else
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                            temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p(:,1,tP,k);zeros(Nb,1)],'l'),aj_p(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_pk,w_ph|y1);
                        cov_wpk_wph_y1=cell(K,K);
                        for h=1:K
                            for k=h+1:K
                                cov_wpk_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(cov_wpk_y1z,1);
                                    if not(blocks(end)==size(cov_wpk_y1z,1))
                                        blocks=cat(2,blocks,size(cov_wpk_y1z,1));
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wp_z_y1))
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1))+cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{h}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wp_z_y1))
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h}+cov_wp_z_y1(:,:,t,k)*temp_p{h}(end-rr+1:end,:));
                                    else
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h});
                                    end
                                end
                                temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p(:,1,tP,h);zeros(Nb,1)],'l'),aj_p(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    else
                        cov_wpk_wph_y1=[];
                    end
                else
                    E_wp_y1=[];
                    diag_Var_wp_y1=[];
                    sum_Var_wp_y1=[];
                    M_cov_wb_wp_y1=[];
                    cov_wp_z_y1=[];
                    cov_wpk_wph_y1=[];
                end
                %delete the variables the dimension of which changes every t
                clear temp_p
                clear temp
                if obj.stem_model.stem_data.X_tv
                    sigma_geo=[];
                end
            end
            
            if (obj.stem_model.stem_par.model_type==2||obj.stem_model.stem_par.model_type==3)
                obj.stem_model.stem_data.stem_varset_p.Y{1}=obj.stem_model.stem_data.Y;
            end
            
            ct2=clock;
            disp(['    Conditional E, Var, Cov evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            ct2_estep=clock;
            disp(['  E step ended in ',stem_misc.decode_time(etime(ct2_estep,ct1_estep))]);
            disp('');
        end
        
        function model_changed = M_step(obj,E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,iteration)
            %DESCRIPTION: M-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1)
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            %sigma_eps                      - [double]          (NxN) sigma_eps
            %st_kalmansmoother_result       - [st_kalmansmoother_result object] (1x1)
            %iteration                      - [double]          (1x1) EM iteration number
            %model_changed                  - [boolean]         (1x1) 1: the number of model parameters changed from the previous iteration (can happen with clustering); 0: no change
            %
            %OUTPUT
            %none: the stem_par property of the stem_model object is updated
            
            disp('  M step started...');
            ct1_mstep=clock;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            K=obj.stem_model.stem_par.k;
            M=obj.stem_model.stem_data.M;
            dim=obj.stem_model.stem_data.dim;
            
            par=obj.stem_model.stem_par;
            st_par_em_step=par;
            
            [aj_bp,aj_p,aj_z]=obj.stem_model.get_aj();
            
            d=1./diag(sigma_eps);
            I=1:length(d);
            inv_sigma_eps=sparse(I,I,d);
            
            if (obj.stem_model.stem_par.model_type==2||obj.stem_model.stem_par.model_type==3)
                obj.stem_model.stem_data.stem_varset_p.Y{1}=[];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             beta update                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_beta))
                ct1=clock;
                disp('    beta update started...');
                temp1=zeros(size(obj.stem_model.stem_data.X_beta,2));
                temp2=zeros(size(obj.stem_model.stem_data.X_beta,2),1);
                d=diag(inv_sigma_eps);
                for t=1:T
                    Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                    if obj.stem_model.stem_data.X_beta_tv
                        tB=t;
                    else
                        tB=1;
                    end
                    if size(obj.stem_model.stem_data.X_beta(:,:,tB),1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,tB);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,tB),1),size(obj.stem_model.stem_data.X_beta(:,:,tB),2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,tB);
                    end
                    temp1=temp1+X_beta_orlated(Lt,:)'*stem_misc.D_apply(X_beta_orlated(Lt,:),d(Lt),'l');
                    temp2=temp2+X_beta_orlated(Lt,:)'*stem_misc.D_apply(E_e_y1(Lt,t)+X_beta_orlated(Lt,:)*par.beta,d(Lt),'l');
                end
                st_par_em_step.beta=temp1\temp2;
                ct2=clock;
                disp(['    beta update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              sigma_eps                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('    sigma_eps update started...');
            ct1=clock;
            temp=zeros(N,1);
            temp1=zeros(N,1);
            d=diag(sigma_eps);
            for t=1:T
                Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                %the next two lines are ok only for sigma_eps diagonal
                temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                temp1(~Lt)=d(~Lt);
                temp=temp+temp1;
            end
            temp=temp/T;
            blocks=[0 cumsum(dim)];
            for i=1:length(dim)
                st_par_em_step.sigma_eps(i,i)=mean(temp(blocks(i)+1:blocks(i+1)));
            end
            ct2=clock;
            disp(['    sigma_eps update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %    G and sigma_eta    %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if par.p>0
                if not(obj.stem_model.stem_data.model_type==1)
                    disp('    G and sigma_eta update started...');
                    ct1=clock;
                    if not(obj.stem_model.stem_par.time_diagonal)
                        S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                        S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                        S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                    else
                        S11=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                        S00=diag(diag(st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                        S10=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3)));
                    end
                    
                    temp=S10/S00;
                    if max(eig(temp))<1
                        st_par_em_step.G=temp;
                    else
                        warning('G is not stable. The last G is retained.');
                    end
                    
                    temp=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                    %st_par_em_step.sigma_eta=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                    %st_par_em_step.sigma_eta=(S11-st_par_em_step.G*S10')/T;
                    if min(eig(temp))>0
                        st_par_em_step.sigma_eta=temp;
                    else
                        warning('Sigma eta is not s.d.p. The last s.d.p. solution is retained');
                    end
                    ct2=clock;
                    disp(['    G and sigma_eta update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                else
                    disp('    G update started...');
                    ct1=clock;
                    S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                    
                    temp=zeros(size(par.G));
                    for i=1:par.p
                        temp(i,i)=trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S10))/trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S00));
                    end
                    if max(eig(temp))<1
                        st_par_em_step.G=temp;
                    else
                        warning('G is not stable. The last G is retained.');
                    end
                    ct2=clock;
                    disp(['    G update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                    
                    disp('    alpha_z update started...');
                    alpha_z=zeros(size(st_par_em_step.alpha_p));
                    for r=1:par.p
                        [~,j_z] = obj.stem_model.get_jz(r);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if obj.stem_model.stem_data.X_bp_tv
                                tBP=t;
                            else
                                tBP=1;
                            end
                            if obj.stem_model.stem_data.X_z_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            
                            if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                            else
                                X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                            end
                            %X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                            
                            %note that here X_z_orlated is not pre-multiplied by aj_z
                            temp1=E_e_y1(:,t)+stem_misc.D_apply(X_z_orlated,aj_z,'l')*stem_misc.D_apply(st_kalmansmoother_result.zk_s(:,t+1),j_z,'l');
                            temp2=(X_z_orlated*stem_misc.D_apply(st_kalmansmoother_result.zk_s(:,t+1),j_z,'l'))';

                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if obj.stem_model.stem_data.model_subtype==1
                                temp1=X_z_orlated;
                                j_z_l=logical(j_z);
                                j_z_l_inv=not(j_z_l);
                                temp1(:,j_z_l_inv)=0;
                                temp2=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                                temp2(:,j_z_l)=0;
                                sum_num=sum_num-trace(temp1*st_kalmansmoother_result.Pk_s(:,:,t+1)*temp2');
                            else
                                %nothing to do as the trace above is zero in this case
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.X_bp))
                                temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                                temp3=zeros(size(temp1,1),1);
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=cat(2,blocks,size(temp1,1));
                                    end
                                    for i=1:length(blocks)-1
                                        temp3(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*temp2(blocks(i)+1:blocks(i+1),:)');
                                    end
                                else
                                    temp3=diag(temp1*temp2');
                                end
                                sum_num=sum_num-sum(temp3(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                    if N>obj.stem_model.system_size
                                        blocks=0:80:size(temp1,1);
                                        if not(blocks(end)==size(temp1,1))
                                            blocks=cat(2,blocks,size(temp1,1));
                                        end
                                        for i=1:length(blocks)-1
                                            temp3(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*temp2(blocks(i)+1:blocks(i+1),:)');
                                        end
                                    else
                                        temp3=diag(temp1*temp2');
                                    end
                                    sum_num=sum_num-sum(temp3(Lt));
                                end
                            end

                            temp1=st_kalmansmoother_result.zk_s(:,t+1)*st_kalmansmoother_result.zk_s(:,t+1)'+st_kalmansmoother_result.Pk_s(:,:,t+1);
                            temp=X_z_orlated;
                            j_z_l=logical(j_z);
                            j_z_l_inv=not(j_z_l);
                            temp(:,j_z_l_inv)=0;
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp,1);
                                if not(blocks(end)==size(temp,1))
                                    blocks=cat(2,blocks,size(temp,1));
                                end
                                for i=1:length(blocks)-1
                                    temp3(blocks(i)+1:blocks(i+1),1)=diag(temp(blocks(i)+1:blocks(i+1),:)*temp1*temp(blocks(i)+1:blocks(i+1),:)');
                                end
                            else
                                temp3=diag(temp*temp1*temp');
                            end
                            sum_den=sum_den+sum(temp3(Lt));
                        end
                        alpha_z(r)=sum_num/sum_den;
                    end
                    st_par_em_step.alpha_z=alpha_z;
                    ct2=clock;
                    disp(['    alpha_z update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                    
                    disp('    v_z and theta_z update started...');
                    ct1=clock;
                    
                    dim=obj.stem_model.stem_data.dim;
                    if obj.stem_model.stem_data.model_subtype==1
                        G_tilde_diag=kron(diag(par.G),ones(dim(1),1));
                    else
                        G_tilde_diag=[];
                        for i=1:par.p
                            G_tilde_diag=cat(1,G_tilde_diag,par.G(i,i)*ones(dim(i),1));
                        end
                    end
                    G_tilde=sparse(1:length(G_tilde_diag),1:length(G_tilde_diag),G_tilde_diag,length(G_tilde_diag),length(G_tilde_diag));
                    temp=S11-S10*G_tilde'-G_tilde*S10'+G_tilde*S00*G_tilde';                    
                    
                    v_temp=par.v_z;
                    %indices are permutated in order to avoid deadlock
                    kindex=randperm(size(par.v_z,1));
                    for k=kindex
                        hindex=randperm(size(par.v_z,1)-k)+k;
                        for h=hindex
                            initial=par.v_z(k,h);
                            ctv1=clock;
                            if obj.stem_model.stem_data.model_subtype==0
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                dim=obj.stem_model.stem_data.stem_varset_p.dim;
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.stem_data.DistMat_z,...
                                    dim(1:par.p),temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            ctv2=clock;
                            disp(['    v_z(',num2str(h),',',num2str(k),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            v_temp(k,h)=min_result;
                            v_temp(h,k)=min_result;
                        end
                    end
                    
                    if min(eig(v_temp))>=0
                        st_par_em_step.v_z=v_temp;
                    else
                        disp('    v_z is not positive definited. The last matrix is retained.');
                    end
                    
                    initial=par.theta_z;
                    ctv1=clock;
                    if obj.stem_model.stem_data.model_subtype==0
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z,par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                    else
                        dim=obj.stem_model.stem_data.stem_varset_p.dim;
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z,par.correlation_type,obj.stem_model.stem_data.DistMat_z,...
                            dim(1:par.p),temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                    end
                    st_par_em_step.theta_z=exp(min_result);
                    ctv2=clock;
                    disp(['    theta_z update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                    
                    ct2=clock;
                    disp(['    v_z and theta_z update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_bp, theta_b and v_b            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_bp))
                disp('    alpha_bp update started...');
                ct1=clock;
                alpha_bp=zeros(size(st_par_em_step.alpha_bp));
                for r=1:obj.stem_model.stem_data.nvar
                    [aj_bp_b,j_b] = obj.stem_model.get_jbp(r);
                    sum_num=0;
                    sum_den=0;
                    for t=1:T
                        if obj.stem_model.stem_data.X_bp_tv
                            tBP=t;
                        else
                            tBP=1;
                        end
                        if obj.stem_model.stem_data.X_z_tv
                            tT=t;
                        else
                            tT=1;
                        end
                        Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                        temp1=E_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp_b,'l');
                        temp2=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t)',M,'r'),obj.stem_model.stem_data.X_bp(:,1,tBP),'r'),j_b,'r');
                        sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                        
                        if par.p>0
                            if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                            else
                                X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                            end
                            X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                            
                            temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),j_b,'l');
                            temp2=zeros(size(temp1,1),1);
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp1,1);
                                if not(blocks(end)==size(temp1,1))
                                    blocks=cat(2,blocks,size(temp1,1));
                                end
                                for i=1:length(blocks)-1
                                    temp2(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                                end
                            else
                                temp2=diag(temp1*X_z_orlated');
                            end
                            sum_num=sum_num-sum(temp2(Lt));
                        end
                        
                        if par.k>0
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            for k=1:K
                                temp1=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),j_b,'l');
                                temp2=[obj.stem_model.stem_data.X_p(:,1,tP,k);zeros(size(temp1,1)-size(obj.stem_model.stem_data.X_p(:,1,tP,k),1),1)];
                                temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',temp2,'r'),aj_p(:,k),'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                        end
                        
                        temp1=E_wb_y1(:,t).^2+diag_Var_wb_y1(:,t);
                        temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(temp1,M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),j_b,'b');
                        sum_den=sum_den+sum(temp1(Lt));
                    end
                    alpha_bp(r,1)=sum_num/sum_den;
                end
                
                st_par_em_step.alpha_bp=alpha_bp;
                
                ct2=clock;
                disp(['    alpha_bp update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                
                disp('    v_b update started...');
                ct1=clock;
                
                if Nb<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wb_y1));
                    for t=1:T
                        temp=temp+E_wb_y1(:,t)*E_wb_y1(:,t)';
                    end
                    temp=temp+sum_Var_wb_y1;
                end
                
                if par.pixel_correlated
                    %indices are permutated in order to avoid deadlock
                    v_temp=par.v_b;
                    kindex=randperm(size(par.v_b,1));
                    for k=kindex
                        hindex=randperm(size(par.v_b,1)-k)+k;
                        for h=hindex
                            initial=par.v_b(k,h);
                            if Nb<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                    obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                    obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            v_temp(k,h)=min_result;
                            v_temp(h,k)=min_result;
                        end
                    end
                    if min(eig(v_temp))>=0
                        st_par_em_step.v_b=v_temp;
                    else
                        disp('    v_b is not positive definited. The last v_b is retained');    
                    end
                    ct2=clock;
                    disp(['    v_b update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                else
                    %nothing because V in this case is the identity matrix
                end
                
                disp('    theta_b updating started...');
                ct1=clock;
                initial=par.theta_b;
                if par.pixel_correlated
                    if Nb<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                            obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        st_par_em_step.theta_b=exp(min_result);
                    else
                        if obj.stem_model.stem_data.stem_varset_b.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_b=exp(min_result);
                        else
                            s=ceil(Nb/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Nb/s);
                            blocks=0:step:Nb;
                            if not(blocks(end)==Nb)
                                blocks=[blocks Nb];
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                end
                                temp=temp+sum_Var_wb_y1(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_b=exp(mean(min_result));
                        end
                    end
                else
                    if Nb<=obj.stem_EM_options.mstep_system_size
                        blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                        for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),...
                                obj.stem_model.stem_data.stem_varset_b.dim(i),temp(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_b(:,i)=exp(min_result);
                        end
                    else
                        blocks_var=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                        for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                            s=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/obj.stem_EM_options.mstep_system_size);
                            step=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/s);
                            blocks=blocks_var(i):step:blocks_var(i+1);
                            if not(blocks(end)==blocks_var(i+1))
                                blocks=cat(2,blocks,blocks_var(i+1));
                            end
                            min_result=[];
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                end
                                temp=temp+sum_Var_wb_y1(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_b(:,i)=exp(mean(min_result));
                        end
                    end
                end
                ct2=clock;
                disp(['    theta_b update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_p               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_p))
                disp('    alpha_p update started...');
                ct1=clock;
                alpha_p=zeros(size(st_par_em_step.alpha_p));
                for s=1:K
                    for r=1:obj.stem_model.stem_data.stem_varset_p.nvar
                        [aj_p_bs,j_p] = obj.stem_model.get_jp(r,s);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if obj.stem_model.stem_data.X_bp_tv
                                tBP=t;
                            else
                                tBP=1;
                            end
                            if obj.stem_model.stem_data.X_z_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            
                            temp1=E_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,s),obj.stem_model.stem_data.X_p(:,1,tP,s),'l'),aj_p_bs,'l');
                            temp2=stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,s)',obj.stem_model.stem_data.X_p(:,1,tP,s),'r'),j_p,'r');
                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if par.p>0
                               
                                if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                    temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                    temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                    X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                                else
                                    X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                                end
                                X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                                                                
                                temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,s),obj.stem_model.stem_data.X_p(:,1,tP,s),'l'),j_p,'l');
                                temp2=zeros(size(temp1,1),1);
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=cat(2,blocks,size(temp1,1));
                                    end
                                    for i=1:length(blocks)-1
                                        temp2(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                                    end
                                else
                                    temp2=diag(temp1*X_z_orlated');
                                end
                                sum_num=sum_num-sum(temp2(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    if not(k==s)
                                        if k<s
                                            kk=s;
                                            ss=k;
                                        else
                                            kk=k;
                                            ss=s;
                                        end
                                        temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{kk,ss}(:,t),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                        temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',[obj.stem_model.stem_data.X_p(:,1,tP,s);zeros(Nb,1)],'r'),j_p,'r');
                                        sum_num=sum_num-sum(temp1(Lt));
                                    end
                                end
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.X_bp))
                                temp1=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,s),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                                temp2=[obj.stem_model.stem_data.X_p(:,1,tP,s);zeros(size(temp1,1)-size(obj.stem_model.stem_data.X_p(:,1,tP,s),1),1)];
                                temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',temp2,'r'),j_p,'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                            
                            temp1=E_wp_y1(:,t,s).^2+diag_Var_wp_y1(:,t,s);
                            temp1=stem_misc.D_apply(stem_misc.D_apply(temp1,obj.stem_model.stem_data.X_p(:,1,tP,s),'b'),j_p,'b');
                            sum_den=sum_den+sum(temp1(Lt));
                        end
                        alpha_p(r,s)=sum_num/sum_den;
                    end
                end
                st_par_em_step.alpha_p=alpha_p;
                ct2=clock;
                disp(['    alpha_p update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                  
                disp('    v_p and theta_p update started...');
                ct1=clock;
                v_temp=par.v_p;
                for z=1:K                    
                    if Np<=obj.stem_EM_options.mstep_system_size
                        temp=zeros(size(sum_Var_wp_y1{z}));
                        for t=1:T
                            temp=temp+E_wp_y1(:,t,z)*E_wp_y1(:,t,z)';
                        end
                        temp=temp+sum_Var_wp_y1{z};
                    end
                    
                    %indices are permutated in order to avoid deadlock
                    kindex=randperm(size(par.v_p(:,:,z),1));
                    for k=kindex
                        hindex=randperm(size(par.v_p(:,:,z),1)-k)+k;
                        for h=hindex
                            initial=par.v_p(k,h,z);
                            ctv1=clock;
                            if Np<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_p(:,:,z),par.theta_p(:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_p(:,:,z),par.theta_p(:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            ctv2=clock;
                            disp(['    v_p(',num2str(h),',',num2str(k),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            v_temp(k,h,z)=min_result;
                            v_temp(h,k,z)=min_result;
                        end
                    end
                    
                    initial=par.theta_p(:,z);
                    ctv1=clock;
                    if Np<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        st_par_em_step.theta_p(:,z)=exp(min_result);
                    else
                        if obj.stem_model.stem_data.stem_varset_p.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_p(:,z)=exp(min_result);
                        else
                            s=ceil(Np/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Np/s);
                            blocks=0:step:Np;
                            if not(blocks(end)==Np)
                                blocks=cat(2,blocks,Np);
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wp_y1(idx,t,z)*E_wp_y1(idx,t,z)';
                                end
                                temp=temp+sum_Var_wp_y1{z}(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_p(:,z)=exp(mean(min_result));
                        end
                        ctv2=clock;
                        disp(['    theta_p(',num2str(z),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                    end
                    ct2=clock;
                    disp(['    v_p and theta_p update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
                for z=1:K
                    if min(eig(v_temp(:,:,z)))>=0
                        st_par_em_step.v_p(:,:,z)=v_temp(:,:,z);
                    else
                        disp(['    v_p(:,:,',num2str(z),') is not positive definited. The last matrix is retained.']);
                    end
                end
            end
            
            model_changed=0;
            if obj.stem_model.stem_par.model_type>=2
                if not(isempty(st_kalmansmoother_result))
                    clear E_e_y1;
                    clear diag_Var_e_y1;
                    clear sigma_eps;
                    clear inv_sigma_eps;
                    clear temp1;
                    clear d;
                    clear I;
                    clear K;
                    if size(obj.stem_model.stem_data.X_z,3)==1
                        %correlation computation
                        for i=1:N
                            L=not(isnan(obj.stem_model.stem_data.Y(i,:)));
                            a=obj.stem_model.stem_data.Y(i,L)';
                            b=st_kalmansmoother_result.zk_s(:,2:end)';
                            b=b(L,:);
                            if not(isempty(b))
                                temp=corr(a,b);
                            else
                                temp=repmat(0.0001,1,par.p);
                            end
                            if obj.stem_model.stem_par.model_type==2
                                temp(temp<=0)=0.0001;
                            end
                            temp(isnan(temp))=0.0001;
                            obj.stem_model.stem_data.X_z(i,:)=temp;
                        end

                        %weight computation
                        for h=1:iteration
                            if obj.stem_model.stem_par.model_type==2
                                obj.stem_model.stem_data.X_z=obj.stem_model.stem_data.X_z.^2;
                                ss=sum(obj.stem_model.stem_data.X_z,2);
                            else
                                obj.stem_model.stem_data.X_z=obj.stem_model.stem_data.X_z.^2.*sign(obj.stem_model.stem_data.X_z);
                                ss=sum(abs(obj.stem_model.stem_data.X_z),2);
                            end
                            for j=1:size(obj.stem_model.stem_data.X_z,2)
                                obj.stem_model.stem_data.X_z(:,j)=obj.stem_model.stem_data.X_z(:,j)./ss;
                            end
                        end
                        
                        %check for empty columns
                        empty=find(sum(abs(obj.stem_model.stem_data.X_z>0.01))==0);
                        if not(isempty(empty))
                            st_par_em_step.p=st_par_em_step.p-length(empty);
                            obj.stem_model.stem_data.X_z(:,empty)=[];
                            sigma_eta_temp=st_par_em_step.sigma_eta;
                            sigma_eta_temp(:,empty)=[];
                            sigma_eta_temp(empty,:)=[];
                            st_par_em_step.sigma_eta=sigma_eta_temp;
                            G_temp=st_par_em_step.G;
                            G_temp(:,empty)=[];
                            G_temp(empty,:)=[]; 
                            st_par_em_step.G=G_temp;
                            model_changed=1;
                            disp(['  Removed ',num2str(length(empty)),' empty cluster(s)']);
                        end
                    else
                            error('X_z must be a 2D matrix and not a 3D array');
                    end
                else
                    error('The Kalman smoother output is empty');
                end
                
                if (obj.stem_model.stem_par.model_type==2||obj.stem_model.stem_par.model_type==3)
                    obj.stem_model.stem_data.stem_varset_p.Y{1}=obj.stem_model.stem_data.Y;
                end
            end
            
            obj.stem_model.stem_par=st_par_em_step;
            ct2_mstep=clock;
            disp(['  M step ended in ',stem_misc.decode_time(etime(ct2_mstep,ct1_mstep))]);
        end
        
        function [E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1] = E_step_parallel(obj,time_steps,st_kalmansmoother_result)
            %DESCRIPTION: parallel version of the E-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1)
            %time_steps                     - [double]          (dTx1) The E-step is computed only for the data related to the time steps in the time_steps vector
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)
            %
            %OUTPUT
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            T=length(time_steps);
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            par=obj.stem_model.stem_par;
            
            fts=time_steps(1);
            
            disp('  E step started...');
            ct1_estep=clock;
            
            [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,~,~,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();

            if p>0
                rr=size(sigma_Z,1);
                if not(obj.stem_model.stem_data.X_z_tv)
                    if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                        temp=obj.stem_model.stem_data.X_z(:,:,1);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,1);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,1),1),size(obj.stem_model.stem_data.X_z(:,:,1),2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    
                    if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                        if obj.stem_model.tapering
                            %is it possible to improve the sparse matrix var_Zt?
                            var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                        else
                            var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                        end
                    else
                        var_Zt=[];
                    end
                end
                if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                    var_Yt=sigma_geo+var_Zt;
                end
            else
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[],[],[]);
                var_Zt=[];
                rr=0;
                %variance of Y
                if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
            end
            
            E_e_y1=obj.stem_model.stem_data.Y(:,time_steps);
            E_e_y1(isnan(E_e_y1))=0;
            if not(isempty(obj.stem_model.stem_data.X_beta))
                disp('    Xbeta evaluation started...');
                ct1=clock;
                Xbeta=zeros(N,T);
                if obj.stem_model.stem_data.X_beta_tv
                    for t=1:T
                        if size(obj.stem_model.stem_data.X_beta(:,:,t+fts-1),1)<N
                            X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,t+fts-1);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,t+fts-1),1),size(obj.stem_model.stem_data.X_beta(:,:,t+fts-1),2))];
                        else
                            X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,t+fts-1);
                        end
                        Xbeta(:,t)=X_beta_orlated*par.beta;
                    end
                else
                    if size(obj.stem_model.stem_data.X_beta(:,:,1),1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,1);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,1),1),size(obj.stem_model.stem_data.X_beta(:,:,1),2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,1);
                    end
                    Xbeta=repmat(X_beta_orlated*par.beta,1,T);
                end
                ct2=clock;
                disp(['    Xbeta evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                E_e_y1=E_e_y1-Xbeta;
            else
                Xbeta=[];
            end
            diag_Var_e_y1=zeros(N,T);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %sigma_Z=Var(Zt)
            %var_Zt=Var(X_z*Zt*X_z')
            
            disp('    Conditional E, Var, Cov evaluation started...');
            ct1=clock;
            
            if not(isempty(obj.stem_model.stem_data.X_bp))
                if obj.stem_model.tapering
                    Lr=find(sigma_W_b);
                    nnz_b=length(Lr);
                end
                %cov_wb_yz time invariant case
                if not(obj.stem_model.stem_data.X_bp_tv)
                    cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp(:,1,1),'r'),aj_bp,'r');
                end
                E_wb_y1=zeros(Nb,T);
                
                if not(obj.stem_model.tapering)
                    sum_Var_wb_y1=zeros(Nb);
                else
                    sum_Var_wb_y1=spalloc(size(sigma_W_b,1),size(sigma_W_b,2),nnz_b);
                end
                
                diag_Var_wb_y1=zeros(Nb,T);
                cov_wb_z_y1=zeros(Nb,rr,T);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_p))
                if obj.stem_model.tapering
                    Lg=find(sigma_W_p{1});
                    nnz_p=length(Lg);
                end
                %cov_wp_yz time invariant case
                if not(obj.stem_model.stem_data.X_p_tv)
                    cov_wp_y=cell(K,1);
                    for k=1:K
                        cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,1,k),'r'),aj_p(:,k),'r');
                    end
                end
                cov_wpk_wph_y1=cell(K,K);
                for h=1:K
                    for k=h+1:K
                        cov_wpk_wph_y1{k,h}=zeros(Np,T);
                    end
                end
                E_wp_y1=zeros(Np,T,K);
                sum_Var_wp_y1=cell(K,1);
                for k=1:K
                    if not(obj.stem_model.tapering)
                        sum_Var_wp_y1{k}=zeros(Np,Np);
                    else
                        sum_Var_wp_y1{k}=spalloc(size(sigma_W_p{k},1),size(sigma_W_p{k},2),nnz_p);
                    end
                end
                diag_Var_wp_y1=zeros(Np,T,K);
                cov_wp_z_y1=zeros(Np,rr,T,K);
            end
            
            if not(isempty(obj.stem_model.stem_data.X_bp)) && not(isempty(obj.stem_model.stem_data.X_p))
                M_cov_wb_wp_y1=zeros(N,T,K);
            else
                M_cov_wb_wp_y1=[];
            end
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(obj.stem_model.stem_data.Y(:,t+fts-1)));
                
                if obj.stem_model.stem_data.X_bp_tv
                    tBP=t+fts-1;
                else
                    tBP=1;
                end
                if obj.stem_model.stem_data.X_z_tv
                    tT=t+fts-1;
                else
                    tT=1;
                end
                if obj.stem_model.stem_data.X_p_tv
                    tP=t+fts-1;
                else
                    tP=1;
                end
                
                %evaluate var_yt in the time variant case
                if obj.stem_model.stem_data.X_tv
                    if not(isempty(obj.stem_model.stem_data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                    end
                    if not(isempty(obj.stem_model.stem_data.X_p))
                        if isempty(sigma_geo)
                            if obj.stem_model.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for k=1:size(obj.stem_model.stem_data.X_p,4)
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,tP,k),'b'),aj_p(:,k),'b');
                        end
                    end
                    if isempty(sigma_geo)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if p>0
                        if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                            temp=obj.stem_model.stem_data.X_z(:,:,tT);
                            temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                        end
                        X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                        
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            if obj.stem_model.tapering
                                %is it possible to improve the sparse matrix var_Zt?
                                var_Zt=sparse(X_z_orlated)*sparse(sigma_Z)*sparse(X_z_orlated');
                            else
                                var_Zt=X_z_orlated*sigma_Z*X_z_orlated';
                            end
                            if isempty(sigma_geo)
                                var_Yt=var_Zt;
                            else
                                var_Yt=sigma_geo+var_Zt;
                            end
                        else
                            var_Zt=[];
                            var_Yt=[];
                        end
                    else
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            var_Yt=sigma_geo;
                        else
                            var_Yt=[];
                        end
                    end
                end
                
                if p>0
                    if N>obj.stem_model.system_size
                        blocks=0:80:size(diag_Var_e_y1,1);
                        if not(blocks(end)==size(diag_Var_e_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                        end
                        for i=1:length(blocks)-1
                            %update diag(Var(e|y1))
                            diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag(X_z_orlated(blocks(i)+1:blocks(i+1),:)*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                        end
                    else
                        temp=X_z_orlated*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                        diag_Var_e_y1(:,t)=diag(temp*X_z_orlated');
                    end
                    %update E(e|y1)
                    temp=st_kalmansmoother_result.zk_s(:,t+fts-1+1);
                    E_e_y1(:,t)=E_e_y1(:,t)-X_z_orlated*temp;
                end
                
                if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), X_z_orlated(Lt,:)*sigma_Z; sigma_Z*X_z_orlated(Lt,:)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    
                    res=obj.stem_model.stem_data.Y(:,time_steps);
                    if not(isempty(Xbeta))
                        res=res-Xbeta;
                    end
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=stem_misc.chol_solve(chol_H1t,temp2(r));
                    else
                        chol_H1t=chol(H1t);
                        cs=stem_misc.chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end

                if not(isempty(obj.stem_model.stem_data.X_bp))
                    %check if the pixel loadings are time variant
                    if obj.stem_model.stem_data.X_bp_tv
                        %cov_wb_yz time variant case
                        cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),obj.stem_model.stem_data.X_bp(:,1,tBP),'r'),aj_bp,'r');
                    end
                    cov_wb_y1z=[cov_wb_y(:,Lt),zeros(size(cov_wb_y,1),rr)];
                    %compute E(w_b|y1);
                    E_wb_y1(:,t)=cov_wb_y1z*cs;
                    %compute Var(w_b|y1)
                    if obj.stem_model.tapering
                        temp_b(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wb_y1z(:,r)');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    else
                        temp_b=stem_misc.chol_solve(chol_H1t,cov_wb_y1z');
                        Var_wb_y1=sigma_W_b-cov_wb_y1z*temp_b;
                    end
                    
                    if p>0
                        %compute cov(w_b,z|y1)
                        cov_wb_z_y1(:,:,t)=temp_b(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                        Var_wb_y1=Var_wb_y1+cov_wb_z_y1(:,:,t)*temp_b(end-rr+1:end,:);
                        %update diag(Var(e|y1))
                        temp=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                        if N>obj.stem_model.system_size
                            blocks=0:80:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                            end
                        else
                            %faster for N
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                        end
                    else
                        cov_wb_z_y1=[];
                    end
                    %compute diag(Var(w_b|y1))
                    diag_Var_wb_y1(:,t)=diag(Var_wb_y1);
                    %compute sum(Var(w_b|y1))
                    sum_Var_wb_y1=sum_Var_wb_y1+Var_wb_y1;
                    %update E(e|y1)
                    E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                    %update diag(Var(e|y1))
                    diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(diag_Var_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                else
                    E_wb_y1=[];
                    diag_Var_wb_y1=[];
                    sum_Var_wb_y1=[];
                    cov_wb_z_y1=[];
                end
                clear temp_b
                if not(isempty(obj.stem_model.stem_data.X_p))
                    %check if the point loadings are time variant
                    if obj.stem_model.stem_data.X_p_tv
                        %cov_wp_yz time invariant case
                        for k=1:K
                            cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_model.stem_data.X_p(:,1,tP,k),'r'),aj_p(:,k),'r');
                        end
                    end
                    temp_p=cell(K,1);
                    for k=1:K
                        cov_wp_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                        %compute E(w_p_k|y1);
                        E_wp_y1(:,t,k)=cov_wp_y1z*cs;
                        %compute Var(w_p_k|y1)
                        if obj.stem_model.tapering
                            temp_p{k}(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wp_y1z(:,r)');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        else
                            temp_p{k}=stem_misc.chol_solve(chol_H1t,cov_wp_y1z');
                            Var_wp_y1=sigma_W_p{k}-cov_wp_y1z*temp_p{k};
                        end
                        
                        if p>0
                            %compute cov(w_p,z|y1)
                            cov_wp_z_y1(:,:,t,k)=temp_p{k}(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+fts-1+1);
                            Var_wp_y1=Var_wp_y1+cov_wp_z_y1(:,:,t,k)*temp_p{k}(end-rr+1:end,:);
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+diag(2*temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %note 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                            end
                        else
                            cov_wp_z_y1=[];
                        end
                        diag_Var_wp_y1(:,t,k)=diag(Var_wp_y1);
                        sum_Var_wp_y1{k}=sum_Var_wp_y1{k}+Var_wp_y1;
                        %update E(e|y1)
                        E_e_y1(:,t)=E_e_y1(:,t)-stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                        %update diag(Var(e|y1))
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(diag_Var_wp_y1(:,t,k),obj.stem_model.stem_data.X_p(:,:,tP,k),'b'),aj_p(:,k),'b');
                        
                        if not(isempty(obj.stem_model.stem_data.X_bp))
                            %compute M_cov(w_b,w_p|y1) namely M*cov(w_b,w_p|y1)
                            if length(M)>obj.stem_model.system_size
                                for i=1:length(M)
                                    if p>0
                                        M_cov_wb_wp_y1(i,t,k)=-cov_wb_y1z(M(i),:)*temp_p{k}(:,i)+cov_wb_z_y1(M(i),:,t)*temp_p{k}(end-rr+1:end,i); %ha gi� l'stem_misc.M_apply su left!!
                                    else
                                        M_cov_wb_wp_y1(i,t,k)=-cov_wb_y1z(M(i),:)*temp_p{k}(:,i);
                                    end
                                end
                            else
                                if p>0
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M))+cov_wb_z_y1(M,:,t)*temp_p{k}(end-rr+1:end,1:length(M))); %ha gi� l'stem_misc.M_apply su left!!
                                else
                                    M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M)));
                                end
                            end
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                            temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p(:,1,tP,k);zeros(Nb,1)],'l'),aj_p(:,k),'l');
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                        end
                    end
                    
                    if K>1
                        %compute cov(w_pk,w_ph|y1);
                        for h=1:K
                            for k=h+1:K
                                cov_wpk_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(cov_wpk_y1z,1);
                                    if not(blocks(end)==size(cov_wpk_y1z,1))
                                        blocks=cat(2,blocks,size(cov_wpk_y1z,1));
                                    end
                                    for i=1:length(blocks)-1
                                        if not(isempty(cov_wp_z_y1))
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1))+cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{h}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                                        else
                                            cov_wpk_wph_y1{k,h}(blocks(i)+1:blocks(i+1),t)=diag(-cov_wpk_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{h}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if not(isempty(cov_wp_z_y1))
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h}+cov_wp_z_y1(:,:,t,k)*temp_p{h}(end-rr+1:end,:));
                                    else
                                        cov_wpk_wph_y1{k,h}(:,t)=diag(-cov_wpk_y1z*temp_p{h});
                                    end
                                end
                                temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[obj.stem_model.stem_data.X_p(:,1,tP,h);zeros(Nb,1)],'l'),aj_p(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    else
                        cov_wpk_wph_y1=[];
                    end
                else
                    E_wp_y1=[];
                    diag_Var_wp_y1=[];
                    sum_Var_wp_y1=[];
                    M_cov_wb_wp_y1=[];
                    cov_wp_z_y1=[];
                    cov_wpk_wph_y1=[];
                end
                clear temp_p
                if obj.stem_model.stem_data.X_tv
                    sigma_geo=[];
                end
            end
            
            ct2=clock;
            disp(['    Conditional E, Var, Cov evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            ct2_estep=clock;
            disp(['  E step ended in ',stem_misc.decode_time(etime(ct2_estep,ct1_estep))]);
            disp('');
        end
        
        function model_changed = M_step_parallel(obj,E_wb_y1,sum_Var_wb_y1,diag_Var_wb_y1,cov_wb_z_y1,E_wp_y1,sum_Var_wp_y1,diag_Var_wp_y1,cov_wp_z_y1,M_cov_wb_wp_y1,cov_wpk_wph_y1,diag_Var_e_y1,E_e_y1,sigma_eps,st_kalmansmoother_result,index,iteration)
            %DESCRIPTION: parallel version of the M-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1)
            %E_wb_y1                        - [double]          (N_bxT) E[wb|Y(1)] conditional expectation of w_b_t with respect to the observed data Y(1)
            %sum_Var_wb_y1                  - [doulbe]          (N_bxN_b) sum(Var[wb|Y(1)]) sum with respect to time of the conditional variance of w_b_t with respect to the observed data
            %diag_Var_wb_y1                 - [double]          (N_bxT) diagonals of Var[wb|Y(1)]
            %cov_wb_z_y1                    - [double]          (N_bxpxT) cov[wb,z_t|Y(1)]
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp|Y(1)]
            %cov_wp_z_y1                    - [double]          (N_pxpxTxK) cov[wp,z|Y(1)]
            %M_cov_wb_wp_y1                 - [double]          (NxTxK)
            %cov_wpk_wph_y1                 - [double]          {KxK}(N_pxT) cov[wp_k,wp_h|Y(1)] k,h=1,...,K
            %diag_Var_e_y1                  - [double]          (NxT) diagonals of Var[e|Y(1)]
            %E_e_y1                         - [double]          (NxT) E[e|Y(1)]
            %sigma_eps                      - [double]          (NxN) sigma_eps
            %st_kalmansmoother_result       - [st_kalmansmoother_result object] (1x1)
            %index                          - [integer >0]      (dKx1) the subset of indices from 1 to K with respect to which estimate the elements of theta_p and v_p
            %iteration                      - [double]          (1x1) EM iteration number
            %model_changed                  - [boolean]         (1x1) 1: the number of model parameters changed from the previous iteration (can happen with clustering); 0: no change
            %
            %OUTPUT
            %none: the stem_par property of the stem_model object is updated
            
            disp('  M step started...');
            ct1_mstep=clock;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            K=obj.stem_model.stem_par.k;
            M=obj.stem_model.stem_data.M;
            dim=obj.stem_model.stem_data.dim;
            
            par=obj.stem_model.stem_par;
            st_par_em_step=par;
            
            [aj_bp,aj_p,aj_z]=obj.stem_model.get_aj();
            
            d=1./diag(sigma_eps);
            I=1:length(d);
            inv_sigma_eps=sparse(I,I,d);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             beta update                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_beta))
                ct1=clock;
                disp('    beta update started...');
                temp1=zeros(size(obj.stem_model.stem_data.X_beta,2));
                temp2=zeros(size(obj.stem_model.stem_data.X_beta,2),1);
                d=diag(inv_sigma_eps);
                for t=1:T
                    Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                    if obj.stem_model.stem_data.X_beta_tv
                        tB=t;
                    else
                        tB=1;
                    end
                    if size(obj.stem_model.stem_data.X_beta(:,:,tB),1)<N
                        X_beta_orlated=[obj.stem_model.stem_data.X_beta(:,:,tB);zeros(N-size(obj.stem_model.stem_data.X_beta(:,:,tB),1),size(obj.stem_model.stem_data.X_beta(:,:,tB),2))];
                    else
                        X_beta_orlated=obj.stem_model.stem_data.X_beta(:,:,tB);
                    end
                    temp1=temp1+X_beta_orlated(Lt,:)'*stem_misc.D_apply(X_beta_orlated(Lt,:),d(Lt),'l');
                    temp2=temp2+X_beta_orlated(Lt,:)'*stem_misc.D_apply(E_e_y1(Lt,t)+X_beta_orlated(Lt,:)*par.beta,d(Lt),'l');
                end
                st_par_em_step.beta=temp1\temp2;
                ct2=clock;
                disp(['    beta update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              sigma_eps                 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('    sigma_eps update started...');
            ct1=clock;
            temp=zeros(N,1);
            temp1=zeros(N,1);
            d=diag(sigma_eps);
            for t=1:T
                Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                %the next two lines are ok only for sigma_eps diagonal
                temp1(Lt)=E_e_y1(Lt,t).^2+diag_Var_e_y1(Lt,t);
                temp1(~Lt)=d(~Lt);
                temp=temp+temp1;
            end
            temp=temp/T;
            blocks=[0 cumsum(dim)];
            for i=1:length(dim)
                st_par_em_step.sigma_eps(i,i)=mean(temp(blocks(i)+1:blocks(i+1)));
            end
            ct2=clock;
            disp(['    sigma_eps update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %    G and sigma_eta    %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if par.p>0
                if not(obj.stem_model.stem_data.model_type==1)
                    disp('    G and sigma_eta update started...');
                    ct1=clock;
                    if not(obj.stem_model.stem_par.time_diagonal)
                        S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                        S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                        S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                    else
                        S11=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                        S00=diag(diag(st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3)));
                        S10=diag(diag(st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'))+diag(diag(sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3)));
                    end
                    
                    temp=S10/S00;
                    if max(eig(temp))<1
                        st_par_em_step.G=temp;
                    else
                        warning('G is not stable. The last G is retained.');
                    end
                    
                    temp=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                    %st_par_em_step.sigma_eta=(S11-S10*par.G'-par.G*S10'+par.G*S00*par.G')/T;
                    %st_par_em_step.sigma_eta=(S11-st_par_em_step.G*S10')/T;
                    if min(eig(temp))>0
                        st_par_em_step.sigma_eta=temp;
                    else
                        warning('Sigma eta is not s.d.p. The last s.d.p. solution is retained');
                    end
                    ct2=clock;
                    disp(['    G and sigma_eta update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                else
                    disp('    G update started...');
                    ct1=clock;
                    S11=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,2:end)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S00=st_kalmansmoother_result.zk_s(:,1:end-1)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.Pk_s(:,:,2:end),3);
                    S10=st_kalmansmoother_result.zk_s(:,2:end)*st_kalmansmoother_result.zk_s(:,1:end-1)'+sum(st_kalmansmoother_result.PPk_s(:,:,2:end),3);
                    
                    temp=zeros(size(par.G));
                    for i=1:par.p
                        temp(i,i)=trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S10))/trace(stem_misc.get_block(dim(1:par.p),i,dim(1:par.p),i,S00));
                    end
                    if max(eig(temp))<1
                        st_par_em_step.G=temp;
                    else
                        warning('G is not stable. The last G is retained.');
                    end
                    ct2=clock;
                    disp(['    G update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                    
                    disp('    alpha_z update started...');
                    alpha_z=zeros(size(st_par_em_step.alpha_p));
                    for r=1:par.p
                        [~,j_z] = obj.stem_model.get_jz(r);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if obj.stem_model.stem_data.X_bp_tv
                                tBP=t;
                            else
                                tBP=1;
                            end
                            if obj.stem_model.stem_data.X_z_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            
                            if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                            else
                                X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                            end
                            %X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                            
                            %note that here X_z_orlated is not pre-multiplied by aj_z
                            temp1=E_e_y1(:,t)+stem_misc.D_apply(X_z_orlated,aj_z,'l')*stem_misc.D_apply(st_kalmansmoother_result.zk_s(:,t+1),j_z,'l');
                            temp2=(X_z_orlated*stem_misc.D_apply(st_kalmansmoother_result.zk_s(:,t+1),j_z,'l'))';

                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if obj.stem_model.stem_data.model_subtype==1
                                temp1=X_z_orlated;
                                j_z_l=logical(j_z);
                                j_z_l_inv=not(j_z_l);
                                temp1(:,j_z_l_inv)=0;
                                temp2=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                                temp2(:,j_z_l)=0;
                                sum_num=sum_num-trace(temp1*st_kalmansmoother_result.Pk_s(:,:,t+1)*temp2');
                            else
                                %nothing to do as the trace above is zero in this case
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.X_bp))
                                temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                                temp3=zeros(size(temp1,1),1);
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=cat(2,blocks,size(temp1,1));
                                    end
                                    for i=1:length(blocks)-1
                                        temp3(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*temp2(blocks(i)+1:blocks(i+1),:)');
                                    end
                                else
                                    temp3=diag(temp1*temp2');
                                end
                                sum_num=sum_num-sum(temp3(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                    if N>obj.stem_model.system_size
                                        blocks=0:80:size(temp1,1);
                                        if not(blocks(end)==size(temp1,1))
                                            blocks=cat(2,blocks,size(temp1,1));
                                        end
                                        for i=1:length(blocks)-1
                                            temp3(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*temp2(blocks(i)+1:blocks(i+1),:)');
                                        end
                                    else
                                        temp3=diag(temp1*temp2');
                                    end
                                    sum_num=sum_num-sum(temp3(Lt));
                                end
                            end

                            temp1=st_kalmansmoother_result.zk_s(:,t+1)*st_kalmansmoother_result.zk_s(:,t+1)'+st_kalmansmoother_result.Pk_s(:,:,t+1);
                            
                            temp=X_z_orlated;
                            j_z_l=logical(j_z);
                            j_z_l_inv=not(j_z_l);
                            temp(:,j_z_l_inv)=0;
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp,1);
                                if not(blocks(end)==size(temp,1))
                                    blocks=cat(2,blocks,size(temp,1));
                                end
                                for i=1:length(blocks)-1
                                    temp3(blocks(i)+1:blocks(i+1),1)=diag(temp(blocks(i)+1:blocks(i+1),:)*temp1*temp(blocks(i)+1:blocks(i+1),:)');
                                end
                            else
                                temp3=diag(temp*temp1*temp');
                            end
                            sum_den=sum_den+sum(temp3(Lt));
                        end
                        alpha_z(r)=sum_num/sum_den;
                    end
                    st_par_em_step.alpha_z=alpha_z;
                    ct2=clock;
                    disp(['    alpha_z update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                    
                    disp('    v_z and theta_z update started...');
                    ct1=clock;
                    
                    dim=obj.stem_model.stem_data.dim;
                    if obj.stem_model.stem_data.model_subtype==1
                        G_tilde_diag=kron(diag(par.G),ones(dim(1),1));
                    else
                        G_tilde_diag=[];
                        for i=1:par.p
                            G_tilde_diag=cat(1,G_tilde_diag,par.G(i,i)*ones(dim(i),1));
                        end
                    end
                    G_tilde=sparse(1:length(G_tilde_diag),1:length(G_tilde_diag),G_tilde_diag,length(G_tilde_diag),length(G_tilde_diag));
                    temp=S11-S10*G_tilde'-G_tilde*S10'+G_tilde*S00*G_tilde';                    
                    
                    v_temp=par.v_z;
                    %indices are permutated in order to avoid deadlock
                    kindex=randperm(size(par.v_z,1));
                    for k=kindex
                        hindex=randperm(size(par.v_z,1)-k)+k;
                        for h=hindex
                            initial=par.v_z(k,h);
                            ctv1=clock;
                            if obj.stem_model.stem_data.model_subtype==0
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                dim=obj.stem_model.stem_data.stem_varset_p.dim;
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_z,par.theta_z,par.correlation_type,obj.stem_model.stem_data.DistMat_z,...
                                    dim(1:par.p),temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            ctv2=clock;
                            disp(['    v_z(',num2str(h),',',num2str(k),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            v_temp(k,h)=min_result;
                            v_temp(h,k)=min_result;
                        end
                    end
                    
                    if min(eig(v_temp))>=0
                        st_par_em_step.v_z=v_temp;
                    else
                        disp('    v_z is not positive definited. The last matrix is retained.');
                    end
                    
                    initial=par.theta_z;
                    ctv1=clock;
                    if obj.stem_model.stem_data.model_subtype==0
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z,par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                    else
                        dim=obj.stem_model.stem_data.stem_varset_p.dim;
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_z,par.correlation_type,obj.stem_model.stem_data.DistMat_z,...
                            dim(1:par.p),temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                    end
                    st_par_em_step.theta_z=exp(min_result);
                    ctv2=clock;
                    disp(['    theta_z update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                    
                    ct2=clock;
                    disp(['    v_z and theta_z update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_bp, theta_b and v_b            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_bp))
                disp('    alpha_bp update started...');
                ct1=clock;
                alpha_bp=zeros(size(st_par_em_step.alpha_bp));
                for r=1:obj.stem_model.stem_data.nvar
                    [aj_bp_b,j_b] = obj.stem_model.get_jbp(r);
                    sum_num=0;
                    sum_den=0;
                    for t=1:T
                        if obj.stem_model.stem_data.X_bp_tv
                            tBP=t;
                        else
                            tBP=1;
                        end
                        if obj.stem_model.stem_data.X_z_tv
                            tT=t;
                        else
                            tT=1;
                        end
                        Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                        temp1=E_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp_b,'l');
                        temp2=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t)',M,'r'),obj.stem_model.stem_data.X_bp(:,1,tBP),'r'),j_b,'r');
                        sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                        
                        if par.p>0
                           
                            if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                            else
                                X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                            end
                            X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                            
                            temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),j_b,'l');
                            temp2=zeros(size(temp1,1));
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(temp1,1);
                                if not(blocks(end)==size(temp1,1))
                                    blocks=cat(2,blocks,size(temp1,1));
                                end
                                for i=1:length(blocks)-1
                                    temp2(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                                end
                            else
                                temp2=diag(temp1*X_z_orlatated');
                            end
                            sum_num=sum_num-sum(temp2(Lt));
                        end
                        
                        if par.k>0
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            for k=1:K
                                temp1=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),j_b,'l');
                                temp2=[obj.stem_model.stem_data.X_p(:,1,tP,k);zeros(size(temp1,1)-size(obj.stem_model.stem_data.X_p(:,1,tP,k),1),1)];
                                temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',temp2,'r'),aj_p(:,k),'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                        end
                        
                        temp1=E_wb_y1(:,t).^2+diag_Var_wb_y1(:,t);
                        temp1=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(temp1,M,'l'),obj.stem_model.stem_data.X_bp(:,1,tBP),'b'),j_b,'b');
                        sum_den=sum_den+sum(temp1(Lt));
                    end
                    alpha_bp(r)=sum_num/sum_den;
                end
                st_par_em_step.alpha_bp=alpha_bp;
                ct2=clock;
                disp(['    alpha_bp update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                
                disp('    v_b update started...');
                ct1=clock;
                if Nb<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wb_y1));
                    for t=1:T
                        temp=temp+E_wb_y1(:,t)*E_wb_y1(:,t)';
                    end
                    temp=temp+sum_Var_wb_y1;
                end
                
                if par.pixel_correlated
                    %indices are permutated in order to avoid deadlock
                    v_temp=par.v_b;
                    kindex=randperm(size(par.v_b,1));
                    for k=kindex
                        hindex=randperm(size(par.v_b,1)-k)+k;
                        for h=hindex
                            initial=par.v_b(k,h);
                            if Nb<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                    obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_b,par.theta_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                    obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            v_temp(k,h)=min_result;
                            v_temp(h,k)=min_result;
                        end
                    end
                    if min(eig(v_temp))>=0
                        st_par_em_step.v_b=v_temp;
                    else
                        disp('    v_b is not positive definited. The last v_b is retained');    
                    end
                    
                    ct2=clock;
                    disp(['    v_b update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
                
                disp('    theta_b updating started...');
                ct1=clock;
                initial=par.theta_b;
                if par.pixel_correlated
                    if Nb<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                            obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        st_par_em_step.theta_b=exp(min_result);
                    else
                        if obj.stem_model.stem_data.stem_varset_b.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b,...
                                obj.stem_model.stem_data.stem_varset_b.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_b=exp(min_result);
                        else
                            s=ceil(Nb/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Nb/s);
                            blocks=0:step:Nb;
                            if not(blocks(end)==Nb)
                                blocks=[blocks Nb];
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                end
                                temp=temp+sum_Var_wb_y1(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_b=exp(mean(min_result));
                        end
                    end
                else
                    if Nb<=obj.stem_EM_options.mstep_system_size
                        blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                        for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),...
                                obj.stem_model.stem_data.stem_varset_b.dim(i),temp(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)),T,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_b(:,i)=exp(min_result);
                        end
                    else
                        blocks_var=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)];
                        for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                            s=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/obj.stem_EM_options.mstep_system_size);
                            step=ceil(obj.stem_model.stem_data.stem_varset_b.dim(i)/s);
                            blocks=blocks_var(i):step:blocks_var(i+1);
                            if not(blocks(end)==blocks_var(i+1))
                                blocks=cat(2,blocks,blocks_var(i+1));
                            end
                            min_result=[];
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wb_y1(idx,t)*E_wb_y1(idx,t)';
                                end
                                temp=temp+sum_Var_wb_y1(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_b,par.correlation_type,obj.stem_model.stem_data.DistMat_b(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_b.tap),log(initial(i)),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_b(:,i)=exp(mean(min_result));
                        end
                    end
                    ct2=clock;
                    disp(['    theta_b update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          alpha_p               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(isempty(obj.stem_model.stem_data.X_p))
                disp('    alpha_p update started...');
                ct1=clock;
                alpha_p=zeros(size(st_par_em_step.alpha_p));
                for s=1:K
                    for r=1:obj.stem_model.stem_data.stem_varset_p.nvar
                        [aj_p_bs,j_p] = obj.stem_model.get_jp(r,s);
                        sum_num=0;
                        sum_den=0;
                        for t=1:T
                            if obj.stem_model.stem_data.X_bp_tv
                                tBP=t;
                            else
                                tBP=1;
                            end
                            if obj.stem_model.stem_data.X_z_tv
                                tT=t;
                            else
                                tT=1;
                            end
                            if obj.stem_model.stem_data.X_p_tv
                                tP=t;
                            else
                                tP=1;
                            end
                            Lt=not(isnan(obj.stem_model.stem_data.Y(:,t)));
                            
                            temp1=E_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,s),obj.stem_model.stem_data.X_p(:,1,tP,s),'l'),aj_p_bs,'l');
                            temp2=stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,s)',obj.stem_model.stem_data.X_p(:,1,tP,s),'r'),j_p,'r');
                            sum_num=sum_num+sum(temp1(Lt).*temp2(Lt)');
                            
                            if par.p>0

                                if (obj.stem_model.stem_data.model_type==1)&&(obj.stem_model.stem_data.model_subtype==0)
                                    temp=obj.stem_model.stem_data.X_z(:,:,tT);
                                    temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                                    X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                                else
                                    X_z_orlated=[obj.stem_model.stem_data.X_z(:,:,tT);zeros(N-size(obj.stem_model.stem_data.X_z(:,:,tT),1),size(obj.stem_model.stem_data.X_z(:,:,tT),2))];
                                end
                                X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');

                                temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,s),obj.stem_model.stem_data.X_p(:,1,tP,s),'l'),j_p,'l');
                                temp2=zeros(size(temp1,1),1);
                                if N>obj.stem_model.system_size
                                    blocks=0:80:size(temp1,1);
                                    if not(blocks(end)==size(temp1,1))
                                        blocks=cat(2,blocks,size(temp1,1));
                                    end
                                    for i=1:length(blocks)-1
                                        temp2(blocks(i)+1:blocks(i+1),1)=diag(temp1(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)');
                                    end
                                else
                                    temp2=diag(temp1*X_z_orlated');
                                end
                                sum_num=sum_num-sum(temp2(Lt));
                            end
                            
                            if K>1
                                for k=1:K
                                    if not(k==s)
                                        if k<s
                                            kk=s;
                                            ss=k;
                                        else
                                            kk=k;
                                            ss=s;
                                        end
                                        temp1=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{kk,ss}(:,t),obj.stem_model.stem_data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                        temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',[obj.stem_model.stem_data.X_p(:,1,tP,s);zeros(Nb,1)],'r'),j_p,'r');
                                        sum_num=sum_num-sum(temp1(Lt));
                                    end
                                end
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.X_bp))
                                temp1=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,s),obj.stem_model.stem_data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                                temp2=[obj.stem_model.stem_data.X_p(:,1,tP,s);zeros(size(temp1,1)-size(obj.stem_model.stem_data.X_p(:,1,tP,s),1),1)];
                                temp1=stem_misc.D_apply(stem_misc.D_apply(temp1',temp2,'r'),j_p,'r');
                                sum_num=sum_num-sum(temp1(Lt));
                            end
                            
                            temp1=E_wp_y1(:,t,s).^2+diag_Var_wp_y1(:,t,s);
                            temp1=stem_misc.D_apply(stem_misc.D_apply(temp1,obj.stem_model.stem_data.X_p(:,1,tP,s),'b'),j_p,'b');
                            sum_den=sum_den+sum(temp1(Lt));
                        end
                        alpha_p(r,s)=sum_num/sum_den;
                    end
                end
                st_par_em_step.alpha_p=alpha_p;
                ct2=clock;
                disp(['    alpha_p update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
                
                %indices are permutated in order to avoid deadlock
                disp('    v_p and theta_p update started...');
                ct1=clock;
                v_temp=par.v_p;
                for z=index %note that z moves over index and not from 1 to K
                    if Np<=obj.stem_EM_options.mstep_system_size
                        temp=zeros(size(sum_Var_wp_y1{z}));
                        for t=1:T
                            temp=temp+E_wp_y1(:,t,z)*E_wp_y1(:,t,z)';
                        end
                        temp=temp+sum_Var_wp_y1{z};
                    end
                    
                    kindex=randperm(size(par.v_p(:,:,z),1));
                    for k=kindex
                        hindex=randperm(size(par.v_p(:,:,z),1)-k)+k;
                        for h=hindex
                            initial=par.v_p(k,h,z);
                            ctv1=clock;
                            if Np<=obj.stem_EM_options.mstep_system_size
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_p(:,:,z),par.theta_p(:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            else
                                disp('WARNING: this operation will take a long time');
                                min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,par.v_p(:,:,z),par.theta_p(:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                    obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            ctv2=clock;
                            disp(['    v_p(',num2str(h),',',num2str(k),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                            v_temp(k,h,z)=min_result;
                            v_temp(h,k,z)=min_result;
                        end
                    end
                    
                    initial=par.theta_p(:,z);
                    ctv1=clock;
                    if Np<=obj.stem_EM_options.mstep_system_size
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        st_par_em_step.theta_p(:,z)=exp(min_result);
                    else
                        if obj.stem_model.stem_data.stem_varset_p.nvar>1
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            st_par_em_step.theta_p(:,z)=exp(min_result);
                        else
                            s=ceil(Np/obj.stem_EM_options.mstep_system_size);
                            step=ceil(Np/s);
                            blocks=0:step:Np;
                            if not(blocks(end)==Np)
                                blocks=cat(2,blocks,Np);
                            end
                            for j=1:length(blocks)-1
                                block_size=blocks(j+1)-blocks(j);
                                idx=blocks(j)+1:blocks(j+1);
                                temp=zeros(block_size);
                                for t=1:T
                                    temp=temp+E_wp_y1(idx,t,z)*E_wp_y1(idx,t,z)';
                                end
                                temp=temp+sum_Var_wp_y1{z}(idx,idx);
                                min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,par.v_p(:,:,z),par.correlation_type,obj.stem_model.stem_data.DistMat_p(idx,idx),...
                                    length(idx),temp,t,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                            end
                            st_par_em_step.theta_p(:,z)=exp(mean(min_result));
                        end
                    end
                    ctv2=clock;
                    disp(['    theta_p(',num2str(z),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                end
                st_par_em_step.v_p=v_temp;
                ct2=clock;
                disp(['    v_p and theta_p update ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end
            
            if (obj.stem_model.stem_par.model_type==2||obj.stem_model.stem_par.model_type==3)
                if not(isempty(st_kalmansmoother_result))
                    clear E_e_y1
                    clear diag_Var_e_y1
                    clear sigma_eps
                    clear inv_sigma_eps
                    clear temp1
                    clear d
                    clear I
                    clear K
                    if size(obj.stem_model.stem_data.X_z,3)==1
                        %correlation computation
                        for i=1:N
                            L=not(isnan(obj.stem_model.stem_data.Y(i,:)));
                            a=obj.stem_model.stem_data.Y(i,L)';
                            b=st_kalmansmoother_result.zk_s(:,2:end)';
                            b=b(L,:);
                            if not(isempty(b))
                                temp=corr(a,b);
                            else
                                temp=repmat(0.0001,1,par.p);
                            end
                            temp(temp<=0)=0.0001;
                            temp(isnan(temp))=0.0001;
                            obj.stem_model.stem_data.X_z(i,:)=temp;
                        end
                        
                        %weight computation
                        for h=1:iteration
                            obj.stem_model.stem_data.X_z=obj.stem_model.stem_data.X_z.^2;
                            ss=sum(obj.stem_model.stem_data.X_z,2);
                            for j=1:size(obj.stem_model.stem_data.X_z,2)
                                obj.stem_model.stem_data.X_z(:,j)=obj.stem_model.stem_data.X_z(:,j)./ss;
                            end
                        end
                    else
                        error('The Kalman smoother output is empty');
                    end
                end
            end
            
            model_changed=0;
            
            obj.stem_model.stem_par=st_par_em_step;
            ct2_mstep=clock;
            disp(['  M step ended in ',stem_misc.decode_time(etime(ct2_mstep,ct1_mstep))]);
        end
        
        function st_par_em_step = M_step_vg_and_theta(obj,E_wp_y1,sum_Var_wp_y1,index)
            %DESCRIPTION: parallel version of the M-step of the EM algorithm only for the parameters v_p and theta_p
            %
            %INPUT
            %obj                            - [stem_EM object]  (1x1)
            %E_wp_y1                        - [double]          (N_pxTxK) E[wp_k|Y(1)]
            %sum_Var_wp_y1                  - [double]          {k}(N_pxN_p) sum(Var[wp_k|Y(1)])
            %diag_Var_wp_y1                 - [double]          (N_pxTxK) diagonals of Var[wp_k|Y(1)]
            %index                          - [integer >0]      (dKx1) the subset of indices from 1 to K with respect to which estimate the elements of theta_p and v_p
            %
            %OUTPUT
            %none: the stem_par property of the stem_model object is updated
            st_par_em_step=obj.stem_model.stem_par;
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            for z=index
                if Np<=obj.stem_EM_options.mstep_system_size
                    temp=zeros(size(sum_Var_wp_y1{z-index(1)+1}));
                    for t=1:size(E_wp_y1,2)
                        temp=temp+E_wp_y1(:,t,z-index(1)+1)*E_wp_y1(:,t,z-index(1)+1)';
                    end
                    temp=temp+sum_Var_wp_y1{z-index(1)+1};
                end
                kindex=randperm(size(st_par_em_step.v_p(:,:,z),1));
                for k=kindex
                    hindex=randperm(size(st_par_em_step.v_p(:,:,z),1)-k)+k;
                    for h=hindex
                        initial=st_par_em_step.v_p(k,h,z);
                        ctv1=clock;
                        if Np<=obj.stem_EM_options.mstep_system_size
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,st_par_em_step.v_p(:,:,z),st_par_em_step.theta_p(:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        else
                            disp('WARNING: this operation will take a long time');
                            min_result = fminsearch(@(x) stem_EM.geo_coreg_function_velement(x,k,h,st_par_em_step.v_p(:,:,z),st_par_em_step.theta_p(:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                                obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),initial,optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        end
                        ctv2=clock;
                        disp(['    v_p(',num2str(h),',',num2str(k),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
                        st_par_em_step.v_p(k,h,z)=min_result;
                        st_par_em_step.v_p(h,k,z)=min_result;
                    end
                end
                
                initial=st_par_em_step.theta_p(:,z);
                ctv1=clock;
                if Np<=obj.stem_EM_options.mstep_system_size
                    min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                        obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                    st_par_em_step.theta_p(:,z)=exp(min_result);
                else
                    if obj.stem_model.stem_data.stem_varset_p.nvar>1
                        disp('WARNING: this operation will take a long time');
                        min_result = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_p,...
                            obj.stem_model.stem_data.stem_varset_p.dim,temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        st_par_em_step.theta_p(:,z)=exp(min_result);
                    else
                        s=ceil(Np/obj.stem_EM_options.mstep_system_size);
                        step=ceil(Np/s);
                        blocks=0:step:Np;
                        if not(blocks(end)==Np)
                            blocks=cat(2,blocks,Np);
                        end
                        for j=1:length(blocks)-1
                            block_size=blocks(j+1)-blocks(j);
                            idx=blocks(j)+1:blocks(j+1);
                            temp=zeros(block_size);
                            for t=1:size(E_wp_y1,2)
                                temp=temp+E_wp_y1(idx,t,z-index(1)+1)*E_wp_y1(idx,t,z-index(1)+1)';
                            end
                            temp=temp+sum_Var_wp_y1{z-index(1)+1}(idx,idx);
                            min_result(j,:) = fminsearch(@(x) stem_EM.geo_coreg_function_theta(x,st_par_em_step.v_p(:,:,z),st_par_em_step.correlation_type,obj.stem_model.stem_data.DistMat_p(idx,idx),...
                                length(idx),temp,obj.stem_model.stem_data.T,obj.stem_model.stem_data.stem_gridlist_p.tap),log(initial),optimset('MaxIter',20,'TolFun',1,'UseParallel','always'));
                        end
                        st_par_em_step.theta_p(:,z)=exp(mean(min_result));
                    end
                end
                ctv2=clock;
                disp(['    theta_p(',num2str(z),') update ended in ',stem_misc.decode_time(etime(ctv2,ctv1))]);
            end
        end
        
        %Class set function
        function set.stem_model(obj,stem_model)
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('You have to provide an object of class stem_model');
            end
        end
    end
    
    
    methods (Static)
        function f = geo_coreg_function_theta(log_theta,v,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to the theta_b or theta_p parameter
            %
            %INPUT
            %log_theta          - [double]      (1x1) natural logarithm of theta
            %v                  - [double]      (qxq) the v_b of v_q matrix
            %correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2  
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') there w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            theta=exp(log_theta);
            n_var=length(var_dims);
            
            if min(eig(v))>0
                if not(isempty(tapering_par))
                    I=zeros(nnz(DistMat),1);
                    J=zeros(nnz(DistMat),1);
                    elements=zeros(nnz(DistMat),1);
                    idx=0;
                    blocks=[0 cumsum(var_dims)];
                    for j=1:n_var
                        for i=j:n_var
                            B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            corr_result=stem_misc.correlation_function(theta,B,correlation_type);
                            weights=stem_misc.wendland(B,tapering_par); %possibile calcolarli una sola volta???
                            corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
                            l=length(corr_result.I);
                            I(idx+1:idx+l)=corr_result.I+blocks(i);
                            J(idx+1:idx+l)=corr_result.J+blocks(j);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                            if not(i==j)
                                I(idx+1:idx+l)=corr_result.J+blocks(j);
                                J(idx+1:idx+l)=corr_result.I+blocks(i);
                                elements(idx+1:idx+l)=corr_result.correlation;
                                idx=idx+l;
                            end
                            
                        end
                    end
                    sigma_W=sparse(I,J,elements);
                else
                    sigma_W=zeros(sum(var_dims));
                    for j=1:n_var
                        for i=j:n_var
                            [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                            if not(isempty(tapering_par))
                                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*stem_misc.wendland(DistMat(block_i,block_j),tapering_par);
                            end
                            if (i~=j)
                                sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                            end
                        end
                    end
                end
                if not(isempty(tapering_par))
                    r = symamd(sigma_W);
                    c=chol(sigma_W(r,r));
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
                else
                    c=chol(sigma_W);
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
                end
            else
                f=10^10;
            end
        end
        
        function f = geo_coreg_function_velement(v_element,row,col,v,theta,correlation_type,DistMat,var_dims,U,T,tapering_par)
            %DESCRIPTION: log-likelihood evaluation with respect to an extra-diagonal element of v_b or v_p
            %
            %INPUT
            %v_element          - [double]      (1x1) the extra-diagonal element of v_b or v_p
            %row                - [double]      (1x1) the row index of the v_element
            %col                - [double]      (1x1) the column index of the v_element
            %v                  - [double]      (qxq) the full v_b or v_p matrix
            %theta              - [double>0]    (1x1) the value of theta_b or theta_p
            %correlation type   - [string]      (1x1) correlation type   - [string]      (1x1) spatial correlation type. 'exponential': exponential spatial correlation function; 'matern32': Matern spatial correlation function with parameter nu=3/2; 'matern52': Matern spatial correlation function with parameter nu=5/2
            %DistMat            - [double]      (N_p x N_p | N_b x N_b) the distance matrix
            %var_dims           - [double]      (qx1) the number of time series for each variable
            %U                  - [double]      (N_p x N_p | N_b x N_b) sum(Var[w|Y(1)]+E[w|Y(1)]*E[w|Y(1)]') there w is w_p or w_b
            %T                  - [integer >0]  (1x1) number of time steps
            %tapering_par       - [double >0]   (1x1) maximum distance after which the spatial correlation is zero
            %
            %OUTPUT
            %f: the log-likelihood value
            
            n_var=length(var_dims);
            v(row,col)=v_element;
            v(col,row)=v_element;
            
            if min(eig(v))>0
                if not(isempty(tapering_par))
                    sigma_W=DistMat;
                else
                    sigma_W=zeros(sum(var_dims));
                end
                
                if not(isempty(tapering_par))
                    I=zeros(nnz(DistMat),1);
                    J=zeros(nnz(DistMat),1);
                    elements=zeros(nnz(DistMat),1);
                    idx=0;
                    blocks=[0 cumsum(var_dims)];
                    for j=1:n_var
                        for i=j:n_var
                            B = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            corr_result=stem_misc.correlation_function(theta,B,correlation_type);
                            weights=stem_misc.wendland(B,tapering_par); %possibile calcolarli una sola volta???
                            corr_result.correlation=v(i,j)*corr_result.correlation.*weights;
                            l=length(corr_result.I);
                            I(idx+1:idx+l)=corr_result.I+blocks(i);
                            J(idx+1:idx+l)=corr_result.J+blocks(j);
                            elements(idx+1:idx+l)=corr_result.correlation;
                            idx=idx+l;
                            if not(i==j)
                                I(idx+1:idx+l)=corr_result.J+blocks(j);
                                J(idx+1:idx+l)=corr_result.I+blocks(i);
                                elements(idx+1:idx+l)=corr_result.correlation;
                                idx=idx+l;
                            end
                        end
                    end
                    sigma_W=sparse(I,J,elements);
                else
                    for j=1:n_var
                        for i=j:n_var
                            [B,block_i,block_j] = stem_misc.get_block(var_dims,i,var_dims,j,DistMat);
                            sigma_W(block_i,block_j)=v(i,j)*stem_misc.correlation_function(theta,B,correlation_type);
                            if not(isempty(tapering_par))
                                sigma_W(block_i,block_j)=sigma_W(block_i,block_j).*stem_misc.wendland(DistMat(block_i,block_j),tapering_par);
                            end
                            if (i~=j)
                                sigma_W(block_j,block_i)=sigma_W(block_i,block_j)';
                            end
                        end
                    end
                end
                if not(isempty(tapering_par))
                    r = symamd(sigma_W);
                    c=chol(sigma_W(r,r));
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(full(c),U(r,r)));
                else
                    c=chol(sigma_W);
                    f=2*T*sum(log(diag(c)))+trace(stem_misc.chol_solve(c,U));
                end
            else
                f=10^10;
            end
        end
    end
end

