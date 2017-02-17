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

classdef stem_krig < handle
    
    %CONSTANTS
    %NN  - number of kriging sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %T   - number of temporal steps
    
    properties
        stem_model=[];  %[stem_model object] (1x1) the stem_model object from which to extract the information for kriging. The model must be estimated.
        X_all=[];       %[double]            (NNxN_bxT) all the loading coefficients evaluated at the kriging sites and for t=1,...,T
        idx_notnan=[];  %[integer >0]        (dNx1) the indices of the non-masked kriging sites (all the kriging sites if they are not masked)
    end
    
    methods
        
        function obj = stem_krig(stem_model)
            %DESCRIPTION: constructor of the class stem_krig
            %
            %INPUT
            %
            %stem_model      - [stem_model object] (1x1) the stem_model object
            %
            %OUTPUT
            %obj             - [stem_krig object]  (1x1) the stem_krig object      
            
            if nargin<1
                error('All the input arguments must be provided');
            end
            obj.stem_model=stem_model;
        end
        
        function st_krig_result = kriging(obj,variable_name,grid,block_size,mask,X,back_transform,no_varcov,crossval)
            %DESCRIPTION: kriging implementation
            %
            %INPUT
            %
            %obj                - [stem_krig object]  (1x1) the stem_krig object 
            %variable_name      - [string]            (1x1) the name of the variable to krige
            %grid               - [stem_grid object]  (1x1) a stem_grid object with the information on the kriging sites
            %<block_size>       - [integer >=0]       (1x1) (default: 0) the maximum number of sites to krige at once. If 0 all the sites are kriged in one step
            %<mask>             - [1 | NaN]           (NNx1) (default: []) a vector mask the elements of which are 1 for the sites that must be kriged and NaN for the sites that do not need kriging
            %<X>                - [double|string]     (NNxN_bxT | 1x1) (default: []) the loading coefficient evaluated at the kriging sites. See Note 2 below.
            %<back_transform>   - [boolean]           (1x1) (default: 1) 1: the kriged variable is back-transformed; 0: no back-transform
            %<no_varcov>        - [boolean]           (1x1) (default: 0) 1: the variance of the kriged variable is not computed; 0: the variance is computed;
            %<crossval>         - [integer>=0]        (1x1) (default: 0) >1: kriging is performed over the cross-validation sites of the variable given by crossval index. This input is reserved to the EM algorithm and should not be used by the user.
            %
            %OUTPUT
            %st_krig_result     - [stem_krig_result object]  (1x1)              
            
            %NOTE 1
            %When the number of kriging sites is large, kriging can be 
            %implemented block by block with dimension given by block_size. 
            %If block_size=0 all the kriging sites are considered in one block.
          
            %NOTE 2
            %X must be provided when the model is based on loading coefficients.
            %X can be either a structure or a string with the path of a folder. 
            %
            %When X is a structure it must have the following variables:
                 %X.X_all       - [double] (NNxN_bxT) the matrix of all the loading coefficients
                 %X.name        - {N_bx1} a cell array with the names of each column of coefficients
                 %X.date_stamp  - (1x1) a stem_datestamp object with the date stamps of the T time steps
            %    
            %When X is a string it must contain the path of the folder
            %where the coefficients are saved block by block. The blocks
            %must contain only the coefficients related to the non-masked
            %sites.
            
            %NOTE 3
            %The structure of each block must be the following
                 %block.data       - (dNxN_bxT) the loading coefficients
                 %block.lat        - (dNx1) the vector of latitudes of the sites of the block
                 %block.lon        - (dNx1) the vector of longitudes of the sites of the block
                 %block.label      - {N_bx1} the name of the N_b loading coefficients
                 %block.idx        - (1x1) the sequence number of the block
                 %block.size       - (1x1) the number of sites of the block
                 %block.date_stamp - (1xT) the date stamps of the T time steps
            
            if nargin<3
                error('Not enough input arguments');
            end
            
            index_var=obj.stem_model.stem_data.stem_varset_p.get_Y_index(variable_name);
            if isempty(index_var)
                error('The variable name is incorrect');
            end
            
            if not(isa(grid,'stem_grid'))
                error('The grid input argument must be of class stem_grid');
            end
            
            if nargin>3
                if block_size<0||isempty(block_size)
                    error('block_size must be >=0');
                end
            else
                block_size=0;
            end
            
            if nargin<9
                crossval=0;
            end
            
            if nargin>5
                if not(isempty(X))
                    if not(ischar(X)||isstruct(X))
                        error('The input argument X must be either a structure or a string. See note 2 of the class constructor');
                    end
                    if ischar(X)
                        last_char=X(length(X));
                        if not(strcmp(last_char,'/'))
                            X=[X,'/'];
                        end
                    end
                    if isstruct(X)
                        if not(isfield(X,'X_all'))
                            error('The field X_all of struct X is missing');
                        end
                        if not(isfield(X,'name'))
                            error('The field name of struct X is missing');
                        end
                        if not(isfield(X,'date_stamp'))
                            error('The field date_stamp of struct X is missing');
                        end
                    end
                else
                    if crossval==0
                        error('The input argument X must be provided');
                    end
                end
            end
            
            if nargin<5
                mask=[];
                obj.idx_notnan=1:size(grid.coordinate,1);
                obj.idx_notnan=obj.idx_notnan';
            else
                if not(isempty(mask))
                    if not(strcmp(grid.grid_type,'regular'))
                        error('mask can be provided only in the case of a regular grid');
                    end
                    if size(mask,2)>1
                        error('mask must be a vector');
                    end
                    if not(length(mask)==size(grid.coordinate,1))
                        error('mask must be a vector with the dimension of the kriging sites');
                    end
                    obj.idx_notnan=find(not(isnan(mask)));
                else
                    obj.idx_notnan=1:size(grid.coordinate,1);
                    obj.idx_notnan=obj.idx_notnan';
                end
            end
            if nargin<6
                obj.X_all=[];
            end
            if (nargin<7)||isempty(back_transform)
                back_transform=1;
            end
            if nargin<8
                no_varcov=0;
            end

            
            if (crossval>0)&&not(obj.stem_model.cross_validation)
                error('The stem_model object does not contain cross-validation information');
            end
            if (crossval>0)&&(not(isempty(X)))
                disp('WARNING: the X provided is not considered as the covariates of cross validation are used');
            end
            if (crossval>0)&&(not(isempty(grid)))
                disp('WARNING: the grid provides in not considered as the grid of cross validation is used');
            end                
            %             if crossval
            %                 grid=obj.stem_model.stem_data.stem_crossval.stem_gridlist.grid{1};
            %             end
            
            if block_size==0
                blocks_krig=[0 size(obj.idx_notnan,1)];
            else
                blocks_krig=0:block_size:size(obj.idx_notnan,1);
                if not(blocks_krig(end)==size(obj.idx_notnan,1))
                    blocks_krig=[blocks_krig size(obj.idx_notnan,1)];
                end
            end
            
            stem_datestamp=obj.stem_model.stem_data.stem_datestamp;
            
            if crossval==0
                %indexes recovering and block test
                loadfromfile=0;
                if not(isempty(X))
                    if isstruct(X)
                        %the X are directly provided
                        if not(isempty(X.X_all))
                            if not(size(X.X_all,1)==size(grid.coordinate,1))
                                error('The size of X_all does not match the number of coordinates');
                            end
                            obj.X_all=X.X_all(obj.idx_notnan,:,:);
                        end
                        
                        idx_datestamp=[];
                        for j=1:length(stem_datestamp.stamp)
                            idx=find(X.date_stamp.stamp==stem_datestamp.stamp(j),1);
                            if isempty(idx)
                                error('The kriging block does not contain the correct datestamp');
                            end
                            idx_datestamp=cat(2,idx_datestamp,idx);
                        end
                        
                        %recover the indices
                        if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_bp_name))
                            X_bp_name=obj.stem_model.stem_data.stem_varset_p.X_bp_name{index_var};
                            if not(isempty(X_bp_name))
                                idx_bp=[];
                                for j=1:length(X_bp_name)
                                    cmp=strcmp(X.name,X_bp_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_bp=cat(2,idx_bp,find(cmp,1));
                                    end
                                end
                            else
                                idx_bp=[];
                            end
                        else
                            idx_bp=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_p_name))
                            X_p_name=obj.stem_model.stem_data.stem_varset_p.X_p_name{index_var};
                            if not(isempty(X_p_name))
                                idx_p=[];
                                for j=1:length(X_p_name)
                                    cmp=strcmp(X.name,X_p_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_p=cat(2,idx_p,find(cmp,1));
                                    end
                                end
                            else
                                idx_p=[];
                            end
                        else
                            idx_p=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta_name))
                            X_beta_name=obj.stem_model.stem_data.stem_varset_p.X_beta_name{index_var};
                            if not(isempty(X_beta_name))
                                idx_beta=[];
                                for j=1:length(X_beta_name)
                                    cmp=strcmp(X.name,X_beta_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idx_beta=cat(2,idx_beta,find(cmp,1));
                                    end
                                end
                            else
                                idx_beta=[];
                            end
                        else
                            idx_beta=[];
                        end
                        
                        if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z_name))
                            X_z_name=obj.stem_model.stem_data.stem_varset_p.X_z_name{index_var};
                            if not(isempty(X_z_name))
                                idX_z=[];
                                for j=1:length(X_z_name)
                                    cmp=strcmp(X.name,X_z_name{j});
                                    if sum(cmp)==0
                                        error('The block does not include the requested covariates');
                                    else
                                        idX_z=cat(2,idX_z,find(cmp,1));
                                    end
                                end
                            else
                                idX_z=[];
                            end
                        else
                            idX_z=[];
                        end
                    else
                        %a folder name has been provided
                        loadfromfile=1;
                        folder=X;
                        files=dir([X,'*.mat']);
                        if isempty(files)
                            error(['The directory ',X,' does not contain the kriging blocks']);
                        end
                        if block_size==0
                            warning('The block_size input argument should be >0. The block_size is recovered from the first block saved on file.');
                            load([X,files(1).name]);
                            block_size=block.size;
                            disp(['block_size is now equal to: ',num2str(block_size)]);
                            %recomputing blocks_krig
                            blocks_krig=0:block_size:size(obj.idx_notnan,1);
                            if not(blocks_krig(end)==size(obj.idx_notnan,1))
                                blocks_krig=[blocks_krig size(obj.idx_notnan,1)];
                            end
                        end
                        %check blocks
                        disp('Hard disk kriging blocks evaluation started...');
                        counter=0;
                        for i=1:length(files)
                            block=[];
                            load([X,files(i).name]); %load X_krig_block variable
                            idx_datestamp=[];
                            for j=1:length(stem_datestamp.stamp)
                                idx=find(block.date_stamp==stem_datestamp.stamp(j),1);
                                if isempty(idx)
                                    error('The kriging block does not contain the correct datestamp');
                                end
                                idx_datestamp=cat(2,idx_datestamp,idx);
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_bp_name))
                                X_bp_name=obj.stem_model.stem_data.stem_varset_p.X_bp_name{index_var};
                                if not(isempty(X_bp_name))
                                    idx_bp=[];
                                    for j=1:length(X_bp_name)
                                        cmp=strcmp(block.label,X_bp_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_bp=cat(2,idx_bp,find(cmp,1));
                                        end
                                    end
                                else
                                    idx_bp=[];
                                end
                            else
                                idx_bp=[];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_p_name))
                                X_p_name=obj.stem_model.stem_data.stem_varset_p.X_p_name{index_var};
                                if not(isempty(X_p_name))
                                    idx_p=[];
                                    for j=1:length(X_p_name)
                                        cmp=strcmp(block.label,X_p_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_p=cat(2,idx_p,find(cmp,1));
                                        end
                                    end
                                else
                                    idx_p=[];
                                end
                            else
                                idx_p=[];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta_name))
                                X_beta_name=obj.stem_model.stem_data.stem_varset_p.X_beta_name{index_var};
                                if not(isempty(X_beta_name))
                                    idx_beta=[];
                                    for j=1:length(X_beta_name)
                                        cmp=strcmp(block.label,X_beta_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idx_beta=cat(2,idx_beta,find(cmp,1));
                                        end
                                    end
                                else
                                    idx_beta=[];
                                end
                            else
                                idx_beta=[];
                            end
                            
                            if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z_name))
                                X_z_name=obj.stem_model.stem_data.stem_varset_p.X_z_name{index_var};
                                if not(isempty(X_z_name))
                                    idX_z=[];
                                    for j=1:length(X_z_name)
                                        cmp=strcmp(block.label,X_z_name{j});
                                        if sum(cmp)==0
                                            error('The block does not include the requested covariates');
                                        else
                                            idX_z=cat(2,idX_z,find(cmp,1));
                                        end
                                    end
                                else
                                    idX_z=[];
                                end
                            else
                                idX_z=[];
                            end
                            
                            if not(size(block.data,1)==blocks_krig(i+1)-blocks_krig(i))
                                error(['Wrong number of rows in kriging block ',num2str(i),'. The block_size input argument must be equal to the number of rows in each block loaded from file.']);
                            end
                            counter=counter+size(block.data,1);
                        end
                        if not(counter==size(obj.idx_notnan,1))
                            error('The total number of rows in kriging blocks is wrong');
                        end
                        disp('Hard disk kriging blocks evaluation ended');
                    end
                end
            else
                %kriging using cross-validation data
            end
            
            
            disp('Kriging started...');
            K=obj.stem_model.stem_par.k;
            
            st_krig_result=stem_krig_result(variable_name,grid,obj.stem_model.stem_data.stem_gridlist_p.grid{index_var},obj.stem_model.stem_data.shape);
            st_krig_result.y_hat=zeros(size(grid.coordinate,1),obj.stem_model.T);
            if not(no_varcov)
                st_krig_result.diag_Var_y_hat=zeros(size(grid.coordinate,1),obj.stem_model.T);
            end
            if K>0
                st_krig_result.E_wp_y1=zeros(size(grid.coordinate,1),obj.stem_model.T,K);
            end
            
            for i=1:length(blocks_krig)-1
                ct1=clock;
                disp(['Kriging block ',num2str(i),' of ',num2str(length(blocks_krig)-1)]);
                block_krig=(blocks_krig(i)+1):blocks_krig(i+1);
                block_krig_length=length(block_krig);
                blocks=cumsum(obj.stem_model.stem_data.stem_varset_p.dim);
                
                %Y manage
                Y_add=nan(block_krig_length,obj.stem_model.T);
                obj.stem_model.stem_data.stem_varset_p.Y{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.Y{index_var},Y_add);
                clear Y_add
                
                %X manage
                if crossval==0
                    if not(isempty(obj.X_all))||loadfromfile
                        if loadfromfile
                            load([folder,files(i).name]);
                            block.data=double(block.data(:,:,idx_datestamp));
                        else
                            if size(obj.X_all,3)>1
                                block.data=obj.X_all(block_krig,:,idx_datestamp);
                            else
                                block.data=obj.X_all(block_krig,:,1);
                            end
                            block.lat=grid.coordinate(obj.idx_notnan(block_krig),1);
                            block.lon=grid.coordinate(obj.idx_notnan(block_krig),2);
                        end
                        
                        if not(isempty(idx_bp))
                            X_krig_block=block.data(:,idx_bp,:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_bp_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_p.X_bp_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_bp{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            obj.stem_model.stem_data.stem_varset_p.X_bp{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_bp{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idx_p))
                            X_krig_block=block.data(:,idx_p,:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_p_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_p.X_p_stds{index_var}(j);
                                end
                            end
                            temp=[];
                            for k=1:size(X_krig_block,2)
                                temp1=X_krig_block(:,k,:);
                                temp=cat(4,temp,temp1);
                            end
                            X_krig_block=temp;
                            if size(obj.stem_model.stem_data.stem_varset_p.X_p{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1,:);
                            end
                            obj.stem_model.stem_data.stem_varset_p.X_p{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_p{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idx_beta))
                            X_krig_block=block.data(:,idx_beta,:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_beta_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_p.X_beta_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_beta{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            obj.stem_model.stem_data.stem_varset_p.X_beta{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_beta{index_var},X_krig_block);
                        end
                        
                        if not(isempty(idX_z))
                            X_krig_block=block.data(:,idX_z,:);
                            if obj.stem_model.stem_data.stem_varset_p.standardized
                                for j=1:size(X_krig_block,2)
                                    X_krig_block(:,j,:)=(X_krig_block(:,j,:)-obj.stem_model.stem_data.stem_varset_p.X_z_means{index_var}(j))/obj.stem_model.stem_data.stem_varset_p.X_z_stds{index_var}(j);
                                end
                            end
                            if size(obj.stem_model.stem_data.stem_varset_p.X_z{index_var},3)==1
                                X_krig_block=X_krig_block(:,:,1);
                            end
                            obj.stem_model.stem_data.stem_varset_p.X_z{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_z{index_var},X_krig_block);
                        end
                    else
                        block.lat=grid.coordinate(obj.idx_notnan(block_krig),1);
                        block.lon=grid.coordinate(obj.idx_notnan(block_krig),2);
                    end
                else
                    %cross-validation data
                    block.lat=grid.coordinate(obj.idx_notnan(block_krig),1);
                    block.lon=grid.coordinate(obj.idx_notnan(block_krig),2);
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_bp))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_bp{1}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_bp{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_bp{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_p))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_p{1}(obj.idx_notnan(block_krig),:,:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_p{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_p{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_beta))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_beta{1}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_beta{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_beta{index_var},X_krig_block);
                    end
                    
                    if not(isempty(obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_z))
                        X_krig_block=obj.stem_model.stem_data.stem_crossval.stem_varset{crossval}.X_z{1}(obj.idx_notnan(block_krig),:,:);
                        obj.stem_model.stem_data.stem_varset_p.X_z{index_var}=cat(1,obj.stem_model.stem_data.stem_varset_p.X_z{index_var},X_krig_block);
                    end
                end

                %Grid manage
                obj.stem_model.stem_data.stem_gridlist_p.grid{index_var}.coordinate=cat(1,obj.stem_model.stem_data.stem_gridlist_p.grid{index_var}.coordinate,[block.lat,block.lon]);
                obj.stem_model.stem_data.update_data();
                obj.stem_model.stem_data.update_distance('point');
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end
                
                %kriging
                [y_hat,diag_Var_y_hat,E_wp_y1,diag_Var_wp_y1]=obj.E_step(no_varcov);
                st_krig_result.y_hat(obj.idx_notnan(block_krig),:)=y_hat(blocks(index_var)+1:blocks(index_var)+block_krig_length,:);
                if not(no_varcov)
                    st_krig_result.diag_Var_y_hat(obj.idx_notnan(block_krig),:)=diag_Var_y_hat(blocks(index_var)+1:blocks(index_var)+block_krig_length,:);
                    if K>0
                        st_krig_result.diag_Var_wp_y1(obj.idx_notnan(block_krig),:,:)=diag_Var_wp_y1(blocks(index_var)+1:blocks(index_var)+block_krig_length,:,:);
                    end
                end
                if K>0
                    st_krig_result.E_wp_y1(obj.idx_notnan(block_krig),:,:)=E_wp_y1(blocks(index_var)+1:blocks(index_var)+block_krig_length,:,:);
                end
                
                %restore original
                obj.stem_model.stem_data.stem_varset_p.Y{index_var}(end-block_krig_length+1:end,:)=[];
                if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_bp))
                    obj.stem_model.stem_data.stem_varset_p.X_bp{index_var}(end-block_krig_length+1:end,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_p))
                    obj.stem_model.stem_data.stem_varset_p.X_p{index_var}(end-block_krig_length+1:end,:,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_beta))
                    obj.stem_model.stem_data.stem_varset_p.X_beta{index_var}(end-block_krig_length+1:end,:,:)=[];
                end
                if not(isempty(obj.stem_model.stem_data.stem_varset_p.X_z))
                    obj.stem_model.stem_data.stem_varset_p.X_z{index_var}(end-block_krig_length+1:end,:,:)=[];
                end

                obj.stem_model.stem_data.stem_gridlist_p.grid{index_var}.coordinate(end-block_krig_length+1:end,:)=[];
                obj.stem_model.stem_data.update_data(); 
                obj.stem_model.stem_data.update_distance();
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    obj.stem_model.stem_data.update_M();
                end
                ct2=clock;
                disp(['Kriging block ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
            end

            if back_transform&&(obj.stem_model.stem_data.stem_varset_p.standardized||obj.stem_model.stem_data.stem_varset_p.log_transformed)
                disp('Back-transformation...');
                s=obj.stem_model.stem_data.stem_varset_p.Y_stds{index_var};
                m=obj.stem_model.stem_data.stem_varset_p.Y_means{index_var};
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&not(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    st_krig_result.y_hat=st_krig_result.y_hat*s+m;
                    if not(no_varcov)
                        st_krig_result.diag_Var_y_hat=st_krig_result.diag_Var_y_hat*s^2;
                    end
                end
                
                if (obj.stem_model.stem_data.stem_varset_p.standardized)&&(obj.stem_model.stem_data.stem_varset_p.log_transformed)
                    y_hat=st_krig_result.y_hat;
                    var_y_hat=st_krig_result.diag_Var_y_hat;
                    st_krig_result.y_hat=exp(y_hat*s+m+(var_y_hat*s^2)/2);
                    if not(no_varcov)
                        st_krig_result.diag_Var_y_hat=(exp(var_y_hat*s^2)-1).*exp(2*(y_hat*s+m)+(var_y_hat*s^2));
                    end
                end
                disp('Back-transformation ended.');
            end
            
            if strcmp(grid.grid_type,'regular')
                disp('Data reshaping...');
                st_krig_result.y_hat=reshape(st_krig_result.y_hat,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T);
                if not(no_varcov)
                    st_krig_result.diag_Var_y_hat=reshape(st_krig_result.diag_Var_y_hat,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T);
                end
                if K>0
                    st_krig_result.E_wp_y1=reshape(st_krig_result.E_wp_y1,grid.grid_size(1),grid.grid_size(2),obj.stem_model.T,K);
                end
                disp('Date reshaping ended.');
            end
            
            
            if not(isempty(mask))&&strcmp(grid.grid_type,'regular')
                disp('Applying mask...');
                mask=reshape(mask,grid.grid_size(1),grid.grid_size(2));
                for t=1:size(st_krig_result.y_hat,3)
                    st_krig_result.y_hat(:,:,t)=st_krig_result.y_hat(:,:,t).*mask;
                    if not(no_varcov)
                        st_krig_result.diag_Var_y_hat(:,:,t)=st_krig_result.diag_Var_y_hat(:,:,t).*mask;
                    end
                    for k=1:K
                        st_krig_result.E_wp_y1(:,:,t,k)=st_krig_result.E_wp_y1(:,:,t,k).*mask;
                    end
                end
                disp('Mask applied.');
            end
            
            st_krig_result.variable_name=variable_name;
            st_krig_result.stem_datestamp=obj.stem_model.stem_data.stem_datestamp;
        end
        
        function [y_hat,diag_Var_y_hat,E_wp_y1,diag_Var_wp_y1] = E_step(obj,no_varcov)
            %DESCRIPTION: kriging is based on the E-step of the EM algorithm
            %
            %INPUT
            %obj                            - [stem_krig object]  (1x1)
            %no_varcov                      - [boolean]           (1x1) the variance of the kriged variable is not computed; 0: the variance is computed;
            %
            %OUTPUT
            %y_hat                          - [double]            (NNxT) the variable estimated over the kriging sites
            %diag_Var_y_hat                 - [double]            (NNxT) the variance of the variable estimated over the kriging sites
            %E_wp_y1                        - [double]            (NNxTxK) E[wp|Y(1)] over the kriging sites
            %diag_Var_wp_y1                 - [double]            (NNxTxK) diagonals of Var[wp|Y(1)] over the kriging sites
            
            N=obj.stem_model.stem_data.N;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                Nb=obj.stem_model.stem_data.stem_varset_b.N;
            else
                Nb=0;
            end
            Np=obj.stem_model.stem_data.stem_varset_p.N;
            T=obj.stem_model.stem_data.T;
            K=obj.stem_model.stem_par.k;
            p=obj.stem_model.stem_par.p;
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            if p>0
                if (obj.stem_model.stem_data.model_type==1)
                    st_kalman=stem_kalman(obj.stem_model);
                    [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,aj_p,aj_z,M] = st_kalman.smoother(0,0);
                else
                    [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
                    st_kalmansmoother_result=obj.stem_model.stem_EM_result.stem_kalmansmoother_result;
                end
                rr=size(sigma_Z,1);
                if not(obj.stem_model.stem_data.X_tv)
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
                    end
                    if not(isempty(sigma_geo))&&(not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p)))
                        var_Yt=sigma_geo+var_Zt;
                    end
                end
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
                st_kalmansmoother_result=stem_kalmansmoother_result([],[],[],[],[]);
                var_Zt=[];
                if not(obj.stem_model.stem_data.X_tv)
                    var_Yt=sigma_geo; %sigma_geo includes sigma_eps
                end
                rr=0;
            end
            
            res=data.Y;
            res(isnan(res))=0;
            if not(isempty(data.X_beta))
                Xbeta=zeros(N,T);
                if data.X_beta_tv
                    for t=1:T
                        if size(data.X_beta(:,:,t),1)<N
                            X_beta_orlated=[data.X_beta(:,:,t);zeros(N-size(data.X_beta(:,:,t),1),size(data.X_beta(:,:,t),2))];
                        else
                            X_beta_orlated=data.X_beta(:,:,t);
                        end
                        Xbeta(:,t)=X_beta_orlated*par.beta;
                    end
                else
                    if size(data.X_beta(:,:,1),1)<N
                        X_beta_orlated=[data.X_beta(:,:,1);zeros(N-size(data.X_beta(:,:,1),1),size(data.X_beta(:,:,1),2))];
                    else
                        X_beta_orlated=data.X_beta(:,:,1);
                    end
                    Xbeta=repmat(X_beta_orlated*par.beta,1,T);
                end
                res=res-Xbeta;
                y_hat=Xbeta;
            else
                y_hat=zeros(size(res));
            end
            
            if not(no_varcov)
                diag_Var_e_y1=zeros(N,T);
            else
                diag_Var_e_y1=[];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Conditional expectation, conditional variance and conditional covariance evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if not(isempty(data.X_bp))
                %cov_wb_yz time invariant case
                if not(data.X_bp_tv)
                    cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),data.X_bp(:,1,1),'r'),aj_bp,'r');
                end
                E_wb_y1=zeros(Nb,T);
                diag_Var_wb_y1=zeros(Nb,T);
                cov_wb_z_y1=zeros(Nb,rr,T);
            end

            if not(isempty(data.X_p))
                if not(data.X_p_tv)
                    cov_wp_y=cell(K,1);
                    %cov_wp_yz time invariant case
                    for k=1:K
                        cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p(:,1,1,k),'r'),aj_p(:,k),'r');
                    end
                end
                cov_wpk_wph_y1=cell(K,K);
                for h=1:K
                    for k=h+1:K
                        cov_wpk_wph_y1{k,h}=zeros(Np,T);
                    end
                end
                E_wp_y1=zeros(Np,T,K);
                if not(no_varcov)
                    diag_Var_wp_y1=zeros(Np,T,K);
                else
                    diag_Var_wp_y1=[];
                end
                cov_wp_z_y1=zeros(Np,rr,T,K);
            else
                E_wp_y1=[];
                diag_Var_wp_y1=[];
            end
            
            if not(isempty(data.X_bp)) && not(isempty(data.X_p))
                M_cov_wb_wp_y1=zeros(N,T,K);
            else
                M_cov_wb_wp_y1=[];
            end
            
            for t=1:T
                %missing at time t
                Lt=not(isnan(data.Y(:,t)));
                
                if data.X_bp_tv
                    tBP=t;
                else
                    tBP=1;
                end
                if data.X_z_tv
                    tT=t;
                else
                    tT=1;
                end
                if data.X_p_tv
                    tP=t;
                else
                    tP=1;
                end
                
                %evaluate var_yt in the time variant case
                if data.X_tv
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
                            var_Yt=[];
                            var_Zt=[];
                        end
                    else
                        if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                            var_Yt=sigma_geo;
                        else
                            var_Yt=[];
                        end
                    end
                end
                
                %check if the temporal loadings are time variant
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
                    y_hat(:,t)=y_hat(:,t)+X_z_orlated*temp;
                end
                
                if not(isempty(obj.stem_model.stem_data.X_bp))||not(isempty(obj.stem_model.stem_data.X_p))
                    %build the Ht matrix
                    if not(isempty(var_Zt))
                        H1t=[var_Yt(Lt,Lt), X_z_orlated(Lt,:)*sigma_Z; sigma_Z*X_z_orlated(Lt,:)', sigma_Z];
                    else
                        H1t=var_Yt(Lt,Lt);
                        temp=[];
                    end
                    
                    if obj.stem_model.tapering
                        cs=[];
                        r = symamd(H1t);
                        chol_H1t=chol(H1t(r,r));
                        temp2=[res(Lt,t);temp];
                        cs(r,1)=stem_misc.chol_solve(chol_H1t,temp2(r));
                        clear temp2
                    else
                        chol_H1t=chol(H1t);
                        cs=stem_misc.chol_solve(chol_H1t,[res(Lt,t);temp]);
                    end
                end

                if not(isempty(data.X_bp))
                    %check if the pixel loadings are time variant
                    if data.X_bp_tv
                        %cov_wb_yz time variant case
                        cov_wb_y=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'r'),data.X_bp(:,1,tBP),'r'),aj_bp,'r');
                    end
                    cov_wb_y1z=[cov_wb_y(:,Lt),zeros(size(cov_wb_y,1),rr)];
                    
                    %compute E(w_b|y1);
                    E_wb_y1(:,t)=cov_wb_y1z*cs;
                    
                    if not(no_varcov)
                        %compute diag(Var(w_b|y1))
                        if obj.stem_model.tapering
                            temp_b(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wb_y1z(:,r)');
                        else
                            temp_b=stem_misc.chol_solve(chol_H1t,cov_wb_y1z');
                        end
                        
                        blocks=0:80:size(diag_Var_wb_y1,1);
                        if not(blocks(end)==size(diag_Var_wb_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_wb_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)=diag(sigma_W_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1))-cov_wb_y1z(blocks(i)+1:blocks(i+1),:)*temp_b(:,blocks(i)+1:blocks(i+1)));
                        end
                    end
                    
                    if (p>0)&&(not(no_varcov))
                        %compute cov(w_b,z|y1)
                        cov_wb_z_y1(:,:,t)=temp_b(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                        blocks=0:80:size(diag_Var_wb_y1,1);
                        if not(blocks(end)==size(diag_Var_wb_y1,1))
                            blocks=cat(2,blocks,size(diag_Var_wb_y1,1));
                        end
                        for i=1:length(blocks)-1
                            diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_wb_y1(blocks(i)+1:blocks(i+1),t)+diag(cov_wb_z_y1(blocks(i)+1:blocks(i+1),:,t)*temp_b(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                        end
                        clear temp_b
                        %update diag(Var(e|y1))
                        temp=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(cov_wb_z_y1(:,:,t),M,'l'),data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                        if N>obj.stem_model.system_size
                            blocks=0:80:size(diag_Var_e_y1,1);
                            if not(blocks(end)==size(diag_Var_e_y1,1))
                                blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %notare 2*
                            end
                        else
                            %faster for N small
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                        end
                    else
                        cov_wb_z_y1=[];
                        clear temp_b
                    end
                    %update y_hat
                    y_hat(:,t)=y_hat(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(E_wb_y1(:,t),M,'l'),data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                    %update diag(Var(e|y1))
                    if not(no_varcov)
                        diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(diag_Var_wb_y1(:,t),M,'l'),data.X_bp(:,1,tBP),'b'),aj_bp,'b'); %tested
                    end
                end
                
                if not(isempty(data.X_p))
                    %check if the point loadings are time variant
                    if data.X_p_tv
                        %cov_wp_yz time invariant case
                        for k=1:K
                            cov_wp_y{k}=stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},data.X_p(:,1,tP,k),'r'),aj_p(:,k),'r');
                        end
                    end
                    temp_p=cell(K,1);
                    for k=1:K
                        cov_wp_y1z=[cov_wp_y{k}(:,Lt) zeros(size(cov_wp_y{k},1),rr)];
                        %compute E(w_p_k|y1);
                        E_wp_y1(:,t,k)=cov_wp_y1z*cs;
                        
                        if not(no_varcov)
                            %compute diag(Var(w_p_k|y1))
                            if obj.stem_model.tapering
                                temp_p{k}(r,:)=stem_misc.chol_solve(full(chol_H1t),cov_wp_y1z(:,r)');
                            else
                                temp_p{k}=stem_misc.chol_solve(chol_H1t,cov_wp_y1z');
                            end
                            
                            blocks=0:80:size(diag_Var_wp_y1(:,t,k),1);
                            if not(blocks(end)==size(diag_Var_wp_y1(:,t,k),1))
                                blocks=cat(2,blocks,size(diag_Var_wp_y1(:,t,k),1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(sigma_W_p{k}(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1))-cov_wp_y1z(blocks(i)+1:blocks(i+1),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                            end
                        end
                        
                        if (p>0)&&(not(no_varcov))
                            %compute cov(w_p,z|y1)
                            cov_wp_z_y1(:,:,t,k)=temp_p{k}(end-rr+1:end,:)'*st_kalmansmoother_result.Pk_s(:,:,t+1);
                            blocks=0:80:size(diag_Var_wp_y1(:,t,k),1);
                            if not(blocks(end)==size(diag_Var_wp_y1(:,t,k),1))
                                blocks=cat(2,blocks,size(diag_Var_wp_y1(:,t,k),1));
                            end
                            for i=1:length(blocks)-1
                                diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag_Var_wp_y1(blocks(i)+1:blocks(i+1),t,k)+diag(cov_wp_z_y1(blocks(i)+1:blocks(i+1),:,t,k)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1)));
                            end
                            %update diag(Var(e|y1))
                            temp=stem_misc.D_apply(stem_misc.D_apply(cov_wp_z_y1(:,:,t,k),data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                            if N>obj.stem_model.system_size
                                blocks=0:80:size(diag_Var_e_y1,1);
                                if not(blocks(end)==size(diag_Var_e_y1,1))
                                    blocks=cat(2,blocks,size(diag_Var_e_y1,1));
                                end
                                for i=1:length(blocks)-1
                                    diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)=diag_Var_e_y1(blocks(i)+1:blocks(i+1),t)+2*diag(temp(blocks(i)+1:blocks(i+1),:)*X_z_orlated(blocks(i)+1:blocks(i+1),:)'); %notare 2*
                                end
                            else
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*diag(temp*X_z_orlated');
                            end
                        else
                            cov_wp_z_y1=[];
                        end
                        %y_hat
                        y_hat(:,t)=y_hat(:,t)+stem_misc.D_apply(stem_misc.D_apply(E_wp_y1(:,t,k),data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                        
                        if not(no_varcov)
                            %update diag(Var(e|y1))
                            diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+stem_misc.D_apply(stem_misc.D_apply(diag_Var_wp_y1(:,t,k),data.X_p(:,:,tP,k),'b'),aj_p(:,k),'b'); %K varianze

                            if not(isempty(data.X_bp))
                                %compute M_cov(w_b,w_p|y1); cioè M*cov(w_b,w_p|y1) da tenere in considerazione nelle forme chiuse!
                                if length(M)>obj.stem_model.system_size
                                    blocks=0:80:length(M);
                                    if not(blocks(end)==length(M))
                                        blocks=cat(2,blocks,length(M));
                                    end
                                    for i=1:length(blocks)-1
                                        %tested
                                        if p>0
                                            M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1))+cov_wb_z_y1(M(blocks(i)+1:blocks(i+1)),:,t)*temp_p{k}(end-rr+1:end,blocks(i)+1:blocks(i+1))); %ha gia' l'stem_misc.M_apply su left!!
                                        else
                                            M_cov_wb_wp_y1(blocks(i)+1:blocks(i+1),t,k)=diag(-cov_wb_y1z(M(blocks(i)+1:blocks(i+1)),:)*temp_p{k}(:,blocks(i)+1:blocks(i+1)));
                                        end
                                    end
                                else
                                    if p>0
                                        M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M))+cov_wb_z_y1(M,:,t)*temp_p{k}(end-rr+1:end,1:length(M))); %ha già l'stem_misc.M_apply su left!!
                                    else
                                        M_cov_wb_wp_y1(1:length(M),t,k)=diag(-cov_wb_y1z(M,:)*temp_p{k}(:,1:length(M)));
                                    end
                                end
                                %update diag(Var(e|y1)) - tested
                                temp=stem_misc.D_apply(stem_misc.D_apply(M_cov_wb_wp_y1(:,t,k),data.X_bp(:,1,tBP),'l'),aj_bp,'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[data.X_p(:,1,tP,k);zeros(Nb,1)],'l'),aj_p(:,k),'l');
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    end
                    
                    if (K>1)&&(not(no_varcov))
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
                                temp=stem_misc.D_apply(stem_misc.D_apply(cov_wpk_wph_y1{k,h}(:,t),data.X_p(:,1,tP,k),'l'),aj_p(:,k),'l');
                                temp=stem_misc.D_apply(stem_misc.D_apply(temp,[data.X_p(:,1,tP,h);zeros(Nb,1)],'l'),aj_p(:,h),'l');
                                %update diag(Var(e|y1))
                                diag_Var_e_y1(:,t)=diag_Var_e_y1(:,t)+2*temp;
                            end
                        end
                    end
                    clear temp_p
                end
                if data.X_tv
                    sigma_geo=[];
                end
            end
            diag_Var_y_hat=diag_Var_e_y1;
        end
          
        %Class set methods
        function set.stem_model(obj,stem_model)
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('stem_model must be of class stem_model');
            end
        end
        
        function set.X_all(obj,X_all)
            if size(obj.idx_notnan,1)~=size(X_all,1)
                error('The number of rows of X must be equal to the number of non-masked pixel');
            end
            if sum(isnan(X_all(:)))>0
                error('No NaN are allowed in covariates');
            end
            obj.X_all=X_all;
        end
    end
end