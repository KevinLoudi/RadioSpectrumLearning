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

function daemon(path_distributed_computing)

%DESCRIPTION: code to run on each slave node for distributed computing
%
%INPUT
%path_distributed_computing        - [string] (1x1) path_distributed_computingfull or relative path of the folder to use for distributed computing

%create the temp folder if it does not exist
if not(exist([path_distributed_computing,'temp'],'dir'))
    mkdir([path_distributed_computing,'temp']);
end

h=now;
node_code = round((h-floor(h))*1000000);
disp([datestr(now),' - Machine code: ',num2str(node_code)]);
timeout=70000; %seconds
while(1)
    exit=0;
    disp([datestr(now),' - Waiting request from master...']);
    while not(exit)
        exit=exist([path_distributed_computing,'whoishere.mat'],'file');
        pause(0.1);
    end
    read=0;
    while not(read)
        try
            load([path_distributed_computing,'whoishere.mat']);
            read=1;
            disp([datestr(now),' - Request from the master received.']);
        catch
        end
        pause(0.1);
    end
    machine.node_code=node_code;
    machine.IDrequest=whoishere.IDrequest;
    if exist('st_model','var')
        machine.require_stemmodel=0;
    else
        machine.require_stemmodel=1;
    end
    
    save([path_distributed_computing,'temp/machine_',num2str(node_code),'.mat'],'machine');
    pause(0.5);
    movefile([path_distributed_computing,'temp/machine_',num2str(node_code),'.mat'],[path_distributed_computing,'machine_',num2str(node_code),'.mat']);
    disp([datestr(now),' - Answered the master.']);
    
    if machine.require_stemmodel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           st_model         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([datestr(now),' - Waiting for stem_model object...']);
        exit=0;
        ct1=clock;
        while not(exit)
            exit=exist([path_distributed_computing,'st_model_parallel_',num2str(node_code),'.mat'],'file');
            pause(4);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while waiting for stem_model object.');
                exit=1;
            end
        end
        read=0;
        ct1=clock;
        while not(read)
            try
                load([path_distributed_computing,'st_model_parallel_',num2str(node_code),'.mat']);
                disp([datestr(now),' - stem_model object received.']);
                read=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while loading stem_model object.');
                read=1;
            end
        end
        deleted=0;
        ct1=clock;
        while not(deleted)
            try
                delete([path_distributed_computing,'st_model_parallel_',num2str(node_code),'.mat']);
                disp([datestr(now),' - stem_model object file deleted.']);
                deleted=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while deleting stem_model object.');
                deleted=1;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           st_par           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for stem_par object']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([path_distributed_computing,'st_par_parallel_',num2str(node_code),'.mat'],'file');
        pause(4);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while waiting for stem_par object.');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([path_distributed_computing,'st_par_parallel_',num2str(node_code),'.mat']);
            disp([datestr(now),' - stem_par object received.']);
            read=1;
        catch
        end
        pause(0.05);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while loading stem_par object.');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([path_distributed_computing,'st_par_parallel_',num2str(node_code),'.mat']);
            disp([datestr(now),' - stem_par object file deleted.']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while deleting stem_par object file.');
            deleted=1;
        end
    end
    st_model.stem_par=st_par;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Kalman           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if st_model.stem_par.p>0
        disp([datestr(now),' - Waiting for Kalman filter data...']);
        exit=0;
        ct1=clock;
        while not(exit)
            exit=exist([path_distributed_computing,'kalman_parallel_',num2str(node_code),'.mat'],'file');
            pause(4);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while waiting for Kalman filter data.');
                exit=1;
            end
        end
        read=0;
        ct1=clock;
        while not(read)
            try
                load([path_distributed_computing,'kalman_parallel_',num2str(node_code),'.mat']);
                disp([datestr(now),' - Kalman filter data received.']);
                read=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while loading Kalman filter data.');
                read=1;
            end
        end
        deleted=0;
        ct1=clock;
        while not(deleted)
            try
                delete([path_distributed_computing,'kalman_parallel_',num2str(node_code),'.mat']);
                disp([datestr(now),' - Kalman filter data file deleted.']);
                deleted=1;
            catch
            end
            pause(0.1);
            ct2=clock;
            if etime(ct2,ct1)>timeout
                disp('Timeout while deleting Kalman filter data.');
                deleted=1;
            end
        end
        
        st_kalman=stem_kalman(st_model);
        st_kalman.filter(0,0,data.time_steps,path_distributed_computing);
        clear data
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         EM -  E step       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for E-step job...']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([path_distributed_computing,'data_parallel_',num2str(node_code),'.mat'],'file');
        pause(4);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while waiting for E-step job.');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([path_distributed_computing,'data_parallel_',num2str(node_code),'.mat']);
            disp([datestr(now),' - E-step job received.']);
            read=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while loading E-step job data.');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([path_distributed_computing,'data_parallel_',num2str(node_code),'.mat']);
            disp([datestr(now),' - E-step job data file deleted.']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while deleting E-step job data file.');
            deleted=1;
        end
    end
    
    
    st_EM_options=stem_EM_options(0.001,1,'single',[],0,[]);
    st_EM=stem_EM(st_model,st_EM_options);
    ct1=clock;
    [output.E_wb_y1,output.sum_Var_wb_y1,output.diag_Var_wb_y1,output.cov_wb_z_y1,output.E_wp_y1,...
        output.sum_Var_wp_y1,output.diag_Var_wp_y1,output.cov_wp_z_y1,output.M_cov_wb_wp_y1,...
        output.cov_wpk_wph_y1,output.diag_Var_e_y1,output.E_e_y1] = st_EM.E_step_parallel(data.time_steps,data.st_kalmansmoother_result);
    output.sum_Var_wb_y1=triu(output.sum_Var_wb_y1);
    for k=1:length(output.sum_Var_wp_y1)
        output.sum_Var_wp_y1{k}=triu(output.sum_Var_wp_y1{k});
    end
    ct2=clock;
    output.ct=etime(ct2,ct1);
    output.cb=data.cb;
    output.iteration=data.iteration;
    output.time_steps=data.time_steps;
    output.node_code=node_code;
    disp([datestr(now),' - Saving E-step job result...']);
    save([path_distributed_computing,'temp/output_',num2str(node_code),'.mat'],'output','-v7.3');
    pause(0.5);
    movefile([path_distributed_computing,'temp/output_',num2str(node_code),'.mat'],[path_distributed_computing,'output_',num2str(node_code),'.mat']);
    disp([datestr(now),' - E-step job result saved.']);
    clear data
    clear output
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         EM - M-Step        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp([datestr(now),' - Waiting for M-step job...']);
    exit=0;
    ct1=clock;
    while not(exit)
        exit=exist([path_distributed_computing,'data_parallel_mstep',num2str(node_code),'.mat'],'file');
        pause(4);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while waiting for M-step job.');
            exit=1;
        end
    end
    read=0;
    ct1=clock;
    while not(read)
        try
            load([path_distributed_computing,'data_parallel_mstep',num2str(node_code),'.mat']);
            disp([datestr(now),' - M-step job received.']);
            read=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while loading M-step job.');
            read=1;
        end
    end
    deleted=0;
    ct1=clock;
    while not(deleted)
        try
            delete([path_distributed_computing,'data_parallel_mstep',num2str(node_code),'.mat']);
            disp([datestr(now),' - M-step job data file deleted.']);
            deleted=1;
        catch
        end
        pause(0.1);
        ct2=clock;
        if etime(ct2,ct1)>timeout
            disp('Timeout while deleting M-step job data file.');
            deleted=1;
        end
    end
    
    output.iteration=data.iteration;
    output.node_code=node_code;
    if not(isempty(data.index))
        disp([datestr(now),' - M-step running...']);
        ct1=clock;
        output.mstep_par=st_EM.M_step_vg_and_theta(data.E_wp_y1,data.sum_Var_wp_y1,data.index);
        ct2=clock;
        output.ct=etime(ct2,ct1);
        output.index=data.index;
        disp([datestr(now),' - Saving M-step job result...']);
    else
        output.index=[];
        disp([datestr(now),' - Nothing to do. Saving empty M-step job result...']);
    end
    save([path_distributed_computing,'temp/output_mstep_',num2str(node_code),'.mat'],'output');
    pause(0.5);
    movefile([path_distributed_computing,'temp/output_mstep_',num2str(node_code),'.mat'],[path_distributed_computing,'output_mstep_',num2str(node_code),'.mat']);
    disp([datestr(now),' - M-step result saved.']);
    clear data
    clear output
end
end
