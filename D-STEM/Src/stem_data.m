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


classdef stem_data < handle
    
    %CONSTANTS
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %N_r = n1_r+...+np_r - total number of elements of the z latent variable whem model_type=1
    %S   = 2 if both point and pixel data are considered. S = 1 if only point data are considered.
    %T   - number of temporal steps
    %TT  = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant
    %p   - dimension of the latent temporal variable z
    
    properties
        stem_varset_p=[];       %[stem_varset object]   (1x1) stem_varset object for the point variables
        stem_varset_b=[];       %[stem_varset object]   (1x1) stem_varset object for the pixel variables
        stem_gridlist_p=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the point variables        
        stem_gridlist_b=[];     %[stem_gridlist object] (1x1) stem_gridlist object for the pixel variables
        stem_datestamp=[]       %[stem_datestamp object](1x1) stem_datestamp object with information on time steps
        stem_crossval=[];       %[stem_crossval object] (1x1) stem_crossval object with information on crossvalidation
        simulated=0;            %[boolean]              (1x1) 1: the data have been simulated; 0: observed data
        
        shape=[];               %[struct]               (1x1) boundary of the geographic region loaded from a shape file

        model_type=0;           %[integer>=0]           (1x1) 0: type 1 model, 1: type 2 model, 2: clustering model
        pixel_correlated=0;     %[boolean]              (1x1) 1: the pixel data are cross-correlated
        
        X_z=[];                 %[double]               (NxpxTT) if model type=0; (N x 1 x TT) if model_type=1 and model_subtype=0; (N x N_r x TT) if model_type=1 and model_subtype=1. The full X_z matrix
    end
    
    properties (SetAccess = private) 
        Y=[];                   %[double]     (NxT) the full observation matrix
        X_bp=[];                %[double]     (Nx1xTT) the full X_bp matrix
        X_beta=[];              %[double]     (NxN_bxTT) the full X_beta matrix
        X_p=[];                 %[double]     (Nx1xTTxK) the full X_p matrix
        
        DistMat_p=[];           %[double]     (N_pxN_p) distance matrix of the point sites
        DistMat_b=[];           %[double]     (N_bxN_b) distance matrix of the pixel sites
        DistMat_z=[];           %[double]     (N_rxN_r) distance matrix of the latent variable z when model_type=1. It is evaluated only if y and z are not 1:1
        
        M=[];                   %[integer >1] (N_px1) vector of indices of the pixel mapped on the point sites
        
        can_reset=0;            %[boolean]    (1x1) 1: data are saved on disk and can be reloaded; 0: data are only on RAM
    end
    
    properties (Dependent, SetAccess = private)
        X_bp_tv=0;              %[boolean]    (1x1) 1: X_bp is time variant; 0: otherwise
        X_beta_tv=0;            %[boolean]    (1x1) 1: X_beta is time variant; 0: otherwise
        X_z_tv=0;               %[boolean]    (1x1) 1: X_z is time variant; 0: otherwise
        X_p_tv=0;               %[boolean]    (1x1) 1: X_p is time variant; 0: otherwise
        X_tv=0;                 %[boolean]    (1x1) 1: at least one between X_bp, X_beta, X_z and X_p is time variant; 0:otherwise
        
        standardized=0;         %[boolean]    (1x1) 1: Y, X_bp, X_beta, X_z and X_p has been standardized; 0: otherwise
        log_transformed=0;      %[boolean]    (1x1) 1: Y has been log-transformed using the method log_transform; 0: otherwise
        
        model_subtype=[];       %[integer]    (1x1) currently used when model_type=1. 0: X_z{i} has only one column; 1: X_z{i} has more than one column
    end
    
    methods
        
        function obj = stem_data(stem_varset_p,stem_gridlist_p,stem_varset_b,stem_gridlist_b,stem_datestamp,shape,can_reset,stem_crossval,pixel_correlated,model_type)
            %DESCRIPTION: is the constructor of the class stem_data
            %
            %INPUT
            %
            %stem_varset_p      - [stem_varset object]    (1x1) stem_varset object for the point variables
            %stem_gridlist_p    - [stem_gridlist object]  (1x1) stem_gridlist object for the point variables  
            %<stem_varset_b>    - [stem_varset object]    (1x1) (default: []) stem_varset object for the pixel variables
            %<stem_gridlist_b>  - [stem_gridlist object]  (1x1) (default: []) stem_gridlist object for the pixel variables
            %<stem_datestamp>   - [stem_datestamp object] (1x1) (default: []) stem_datestamp object with information on time steps
            %<shape>            - [struct]                      (default: world boundaries) geographic data structure loaded from a shapefile with the boundary of the geographic region
            %<can_reset>        - [boolean]               (1x1) (default: 0) 1: the data are saved on disk and they can be reloaded using the method reset of this class after, for example, data transformation
            %<stem_crossval>    - [stem_crossval object]  (1x1) (default: []) stem_crossval object with information on crossvalidation
            %<model_type>       - [integer>=0]            (1x1) (default: 0) 0: type 1 model, 1: type 2 model, 2: clustering model
            %
            %OUTPUT
            %obj                - [stem_data object]      (1x1) the stem_data object
            
            if nargin<2
                error('Not enough input parameters');
            end
            obj.stem_varset_p=stem_varset_p;
            
            obj.stem_gridlist_p=stem_gridlist_p;
            if not(length(obj.stem_gridlist_p.grid)==length(obj.stem_varset_p.Y))
                error('The number of stem_grids must be equal to the q');
            end
            for i=1:length(obj.stem_gridlist_p.grid)
                if not(size(obj.stem_gridlist_p.grid{i}.coordinate,1)==size(obj.stem_varset_p.Y{i},1))
                    error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                end
                if not(strcmp(obj.stem_gridlist_p.grid{i}.site_type,'point'))
                    error('Only point data are supported in stem_gridlist_p');
                end
            end
            
            if nargin==3
                error('stem_gridlist_b must be provided');
            end
            if nargin>2
                if not(isempty(stem_varset_b))
                    obj.stem_varset_b=stem_varset_b;
                    if not(length(obj.stem_varset_b.dim)==length(stem_varset_p.dim))
                        error('stem_varset_b must contain the same number of variables of stem_varset_p');
                    end
                    if not(size(obj.stem_varset_b.Y{1},2)==size(stem_varset_p.Y{1},2))
                        error('The number of temporal steps in Y cannot differ between stem_varset_b and stem_varset_p');
                    end
                    if not(isempty(obj.stem_varset_b.X_beta))
                        if not(size(obj.stem_varset_b.X_beta{1},3)==size(stem_varset_p.X_beta{1},3))
                            error('The number of temporal steps in X_beta cannot differ between stem_varset_b and stem_varset_p');
                        end
                    end
                    if not(isempty(obj.stem_varset_b.X_z))
                        if not(size(obj.stem_varset_b.X_z{1},3)==size(stem_varset_p.X_z{1},3))
                            error('The number of temporal steps in X_z cannot differ between stem_varset_b and stem_varset_p');
                        end
                    end
                    obj.stem_gridlist_b=stem_gridlist_b;
                    if not(length(obj.stem_gridlist_b.grid)==length(obj.stem_varset_b.Y))
                        error('The number of stem_grids must be equal to the q');
                    end
                    for i=1:length(obj.stem_gridlist_b.grid)
                        if not(size(obj.stem_gridlist_b.grid{i}.coordinate,1)==size(obj.stem_varset_b.Y{i},1))
                            error('The number of coordinates in the grid{i} must be equal to the number of rows of Y{i}');
                        end
                        if not(strcmp(obj.stem_gridlist_b.grid{i}.site_type,'pixel'))
                            error('The grids of stem_gridlist_b must be grids of pixels');
                        end
                        if not(strcmp(obj.stem_gridlist_b.grid{i}.pixel_shape,'square'))
                            error('Only square pixels are supported. Check pixel shape');
                        end
                        if not(obj.stem_gridlist_b.grid{i}.pixel_side_w==obj.stem_gridlist_b.grid{i}.pixel_side_h)
                            error('Only square pixels are supported. Check pixel_side_w and pixel_side_h');
                        end
                    end
                    if not(strcmp(obj.stem_gridlist_b.grid{1}.unit,obj.stem_gridlist_p.grid{1}.unit))
                        error('Both the stem_gridlist objects must contain grids with the same unit');
                    end
                end
            end
            if nargin>4
                obj.stem_datestamp=stem_datestamp;
                if not(length(obj.stem_datestamp.stamp)==obj.stem_varset_p.T)
                    error('The number of datestamps differs from T');
                end
            end

            if nargin>=6
                if not(isempty(shape))
                    obj.shape=shape;
                    obj.shape(1,1).Geometry='Line';
                else
                    obj.shape=shaperead('landareas.shp');
                end
            else
                obj.shape=shaperead('landareas.shp');
            end
            
            if nargin>=7
                if not(isempty(can_reset))
                    obj.can_reset=can_reset;
                else
                    obj.can_reset=0;
                end
            else
                obj.can_reset=0;
            end
            
            if nargin>=8
                if not(isempty(stem_crossval))
                    obj.stem_crossval=stem_crossval;
                    idx_var = obj.stem_varset_p.get_Y_index(obj.stem_crossval.variable_name);
                    if isempty(idx_var)
                        error('Cross-validation variable not found');
                    end
                end
            end
            
            if nargin>=9
                if not(isempty(pixel_correlated))
                    if isempty(obj.stem_varset_b)
                        disp('WARNING: Pixel data are not provided. The pixel_correlated input argument is ignored');
                    else
                        obj.pixel_correlated=pixel_correlated;
                    end
                end
            end
            
            if nargin>=10
                if not(isempty(model_type))
                    obj.model_type=model_type;
                    if obj.model_type==1
                        if not(isempty(obj.stem_varset_p.X_p))
                            error('X_p must be empty when model_type is 1');
                        end
                        if not(isempty(obj.stem_varset_p.X_z))
                            r=size(obj.stem_varset_p.X_z{1},2);
                            for i=2:length(obj.stem_varset_p.X_z)
                                if not(size(obj.stem_varset_p.X_z{i},2)==r)
                                    error('The number of columns of each X_z{i} must be the same when model_type=1');
                                end
                            end
                            n=size(obj.stem_varset_p.X_z{1},1);
                            if size(obj.stem_varset_p.X_z{1},2)>1
                                for i=2:length(obj.stem_varset_p.X_z)
                                    if not(size(obj.stem_varset_p.X_z{i},1)==n)
                                        error('The number of sites of each variable Y{i} must be the same when model_type=1 and X_z{i} has multiple columns');
                                    end
                                end
                            end
                        else
                            error('X_z cannot be empty when model_type=1');
                        end
                    end
                    if (obj.model_type==2||obj.model_type==3)
                        if obj.nvar>1
                            error('The clustering model is only available in the univariate case (q=1)');
                        end
                        if not(isempty(obj.stem_varset_b))
                            error('The clustering model is only available for point level data');
                        end
                        if isempty(obj.stem_varset_p.X_z)
                            error('X_z must be provided when the clustering model is enabled');
                        end
                    end
                end
            end
            
            %if obj.can_reset
            % obj.set_original;
            %end
            
            obj.update_data();
            obj.update_distance();
            if not(isempty(stem_varset_b))
                obj.update_M();
            end
            obj.remove_duplicated_sites();
        end
        
        function update_data(obj)
            %DESCRIPTION: generates the matrices Y, X_bp, X_beta, X_z and X_p 
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            disp('Generating data matrices...');
            %Y
            Y=[];
            for i=1:length(obj.stem_varset_p.Y)
                Y=cat(1,Y,obj.stem_varset_p.Y{i});
            end
            if not(isempty(obj.stem_varset_b))
                for i=1:length(obj.stem_varset_b.Y)
                    Y=cat(1,Y,obj.stem_varset_b.Y{i});
                end
            end
            obj.Y=Y;
            clear Y;
            
            %X_bp
            if not(isempty(obj.stem_varset_b))
                X_bp=[];
                if not(isempty(obj.stem_varset_b.X_bp))&&not(isempty(obj.stem_varset_p.X_bp))
                    done=0;
                    if size(obj.stem_varset_p.X_bp{1},3)==1 && size(obj.stem_varset_b.X_bp{1},3)==obj.T
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,repmat(obj.stem_varset_p.X_bp{i},[1,1,obj.T]));
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_b.X_bp{i});
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_bp{1},3)==obj.T && size(obj.stem_varset_b.X_bp{1},3)==1 && not(done)
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_p.X_bp{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,repmat(obj.stem_varset_b.X_bp{i},[1,1,obj.T]));
                        end
                        done=1;
                    end
                    if size(obj.stem_varset_p.X_bp{1},3)==size(obj.stem_varset_b.X_bp{1},3) && not(done)
                        for i=1:length(obj.stem_varset_p.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_p.X_bp{i});
                        end
                        for i=1:length(obj.stem_varset_b.X_bp)
                            X_bp=cat(1,X_bp,obj.stem_varset_b.X_bp{i});
                        end
                    end
                end
                obj.X_bp=X_bp;
                clear X_bp;
            end

            %X_beta
            if not(isempty(obj.stem_varset_p.X_beta))
                T_max=[];
                for i=1:length(obj.stem_varset_p.X_beta)
                    if not(isempty(obj.stem_varset_p.X_beta{i}))
                        T_max=size(obj.stem_varset_p.X_beta{i},3);
                    end
                end
                if T_max>1
                    n1=0;
                    n2=0;
                    for i=1:length(obj.stem_varset_p.X_beta)
                        n1=n1+size(obj.stem_varset_p.X_beta{i},1);
                        if sum(abs(obj.stem_varset_p.X_beta{i}(:)))>0
                            n2=n2+size(obj.stem_varset_p.X_beta{i},2);
                        end
                    end
                    X_beta=zeros(n1,n2,obj.T);
                    for t=1:size(obj.stem_varset_p.X_beta{1},3)
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_beta)
                            if sum(abs(obj.stem_varset_p.X_beta{i}(:)))>0
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i}(:,:,t));
                            else
                                X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_p.X_beta{i}(:,:,t),1),size(X_temp,2)));
                            end
                        end
                        X_beta(:,:,t)=X_temp;
                    end
                else
                    X_temp=[];
                    for i=1:length(obj.stem_varset_p.X_beta)
                        if sum(abs(obj.stem_varset_p.X_beta{i}(:)))>0
                            X_temp=blkdiag(X_temp,obj.stem_varset_p.X_beta{i});
                        else
                            X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_p.X_beta{i},1),size(X_temp,2)));
                        end
                    end
                    X_beta=X_temp;
                end
                X_beta_p=X_beta;
                clear X_beta;
                clear X_temp;
            else
                X_beta_p=[];
            end
            
            if not(isempty(obj.stem_varset_b))
                if not(isempty(obj.stem_varset_b.X_beta))
                    T_max=[];
                    for i=1:length(obj.stem_varset_b.X_beta)
                        if not(isempty(obj.stem_varset_b.X_beta{i}))
                            T_max=size(obj.stem_varset_b.X_beta{i},3);
                        end
                    end
                    if T_max>1
                        n1=0;
                        n2=0;
                        for i=1:length(obj.stem_varset_b.X_beta)
                            n1=n1+size(obj.stem_varset_b.X_beta{i},1);
                            if sum(abs(obj.stem_varset_b.X_beta{i}(:)))>0
                                n2=n2+size(obj.stem_varset_b.X_beta{i},2);
                            end
                        end
                        X_beta=zeros(n1,n2,obj.T);
                        for t=1:size(obj.stem_varset_b.X_beta{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_b.X_beta)
                                if sum(abs(obj.stem_varset_b.X_beta{i}(:)))>0
                                    X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i}(:,:,t));
                                else
                                    X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_b.X_beta{i}(:,:,t),1),size(X_temp,2)));
                                end
                            end
                            X_beta(:,:,t)=X_temp;
                        end
                    else
                        X_temp=[];
                        for i=1:length(obj.stem_varset_b.X_beta)
                            if sum(abs(obj.stem_varset_b.X_beta{i}(:)))>0
                                X_temp=blkdiag(X_temp,obj.stem_varset_b.X_beta{i});
                            else
                                X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_b.X_beta{i},1),size(X_temp,2)));
                            end
                        end
                        X_beta=X_temp;
                    end
                    X_beta_b=X_beta;
                    clear X_beta;
                    clear X_temp;
                else
                    X_beta_b=[];
                end
            else
                X_beta_b=[];
            end
            
            if not(isempty(X_beta_p))&&not(isempty(X_beta_b))
                obj.X_beta=zeros(size(X_beta_p,1)+size(X_beta_b,1),size(X_beta_p,2)+size(X_beta_b,2),size(X_beta_p,3));
                for t=1:size(X_beta_p,3)
                    obj.X_beta(:,:,t)=blkdiag(X_beta_p(:,:,t),X_beta_b(:,:,t));
                end
            else
              if not(isempty(X_beta_p))
                  obj.X_beta=X_beta_p;
              end
              if not(isempty(X_beta_b))
                  obj.X_beta=cat(1,zeros(obj.stem_varset_p.N,size(X_beta_p,2),obj.T),X_beta_p);
              end
            end
            
            %X_z
            if not(obj.model_type==1)
                if not(isempty(obj.stem_varset_p.X_z))
                    T_max=[];
                    for i=1:length(obj.stem_varset_p.X_z)
                        if not(isempty(obj.stem_varset_p.X_z{i}))
                            T_max=size(obj.stem_varset_p.X_z{i},3);
                        end
                    end
                    if T_max>1
                        n1=0;
                        n2=0;
                        for i=1:length(obj.stem_varset_p.X_z)
                            n1=n1+size(obj.stem_varset_p.X_z{i},1);
                            if sum(abs(obj.stem_varset_p.X_z{i}(:)))>0
                                n2=n2+size(obj.stem_varset_p.X_z{i},2);
                            end
                        end
                        X_z=zeros(n1,n2,obj.T);
                        for t=1:size(obj.stem_varset_p.X_z{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                if sum(abs(obj.stem_varset_p.X_z{i}(:)))>0
                                    X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i}(:,:,t));
                                else
                                    X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_p.X_z{i}(:,:,t),1),size(X_temp,2)));
                                end
                            end
                            X_z(:,:,t)=X_temp;
                        end
                    else
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_z)
                            if sum(abs(obj.stem_varset_p.X_z{i}(:)))>0
                                X_temp=blkdiag(X_temp,obj.stem_varset_p.X_z{i});
                            else
                                X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_p.X_z{i},1),size(X_temp,2)));
                            end
                        end
                        X_z=X_temp;
                    end
                    X_z_p=X_z;
                    clear X_z;
                    clear X_temp;
                else
                    X_z_p=[];
                end
                
                if not(isempty(obj.stem_varset_b))
                    if not(isempty(obj.stem_varset_b.X_z))
                        T_max=[];
                        for i=1:length(obj.stem_varset_b.X_z)
                            if not(isempty(obj.stem_varset_b.X_z{i}))
                                T_max=size(obj.stem_varset_b.X_z{i},3);
                            end
                        end
                        if T_max>1
                            n1=0;
                            n2=0;
                            for i=1:length(obj.stem_varset_b.X_z)
                                n1=n1+size(obj.stem_varset_b.X_z{i},1);
                                if sum(abs(obj.stem_varset_b.X_z{i}(:)))>0
                                    n2=n2+size(obj.stem_varset_b.X_z{i},2);
                                end
                            end
                            X_z=zeros(n1,n2,obj.T);
                            for t=1:size(obj.stem_varset_b.X_z{1},3)
                                X_temp=[];
                                for i=1:length(obj.stem_varset_b.X_z)
                                    if sum(abs(obj.stem_varset_b.X_z{i}(:)))>0
                                        X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i}(:,:,t));
                                    else
                                        X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_b.X_z{i}(:,:,t),1),size(X_temp,2)));
                                    end
                                end
                                X_z(:,:,t)=X_temp;
                            end
                        else
                            X_temp=[];
                            for i=1:length(obj.stem_varset_b.X_z)
                                if sum(abs(obj.stem_varset_b.X_z{i}(:)))>0
                                    X_temp=blkdiag(X_temp,obj.stem_varset_b.X_z{i});
                                else
                                    X_temp=cat(1,X_temp,zeros(size(obj.stem_varset_b.X_z{i},1),size(X_temp,2)));
                                end
                            end
                            X_z=X_temp;
                        end
                        X_z_b=X_z;
                        clear X_z;
                        clear X_temp;
                    else
                        X_z_b=[];
                    end
                else
                    X_z_b=[];
                end
                
                if not(isempty(X_z_p))&&not(isempty(X_z_b))
                    obj.X_z=blkdiag(X_z_p,X_z_b);
                else
                    if not(isempty(X_z_p))
                        obj.X_z=X_z_p;
                    end
                    if not(isempty(X_z_b))
                        obj.X_z=[zeros(obj.stem_varset_p.N,obj.T);X_z_p];
                    end
                end
            else
                %model_type=1
                if size(obj.stem_varset_p.X_z{1},3)==obj.T
                    if obj.model_subtype==1
                        %X_z has more than one column
                        r=obj.dim(1)*size(obj.stem_varset_p.X_z{1},2);
                        X_z=zeros(obj.Np,r,obj.T);
                        for t=1:size(obj.stem_varset_p.X_z{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z)
                                X_temp2=[];
                                for j=1:size(obj.stem_varset_p.X_z{i},2)
                                    X_temp2=cat(2,X_temp2,diag(obj.stem_varset_p.X_z{i}(:,j,t)));
                                end
                                X_temp=cat(1,X_temp,X_temp2);
                            end
                            X_z(:,:,t)=X_temp;
                        end
                    else
                        %X_z has one column
                        X_z=zeros(obj.Np,1,obj.T); %should be a diagonal matrix for each t but it is leaved as a vector to save memory
                        for t=1:size(obj.stem_varset_p.X_z{1},3)
                            X_temp=[];
                            for i=1:length(obj.stem_varset_p.X_z{i})
                                X_temp=cat(1,X_temp,obj.stem_varset_p.X_z{i}(:,t));
                            end
                            X_z(:,1,t)=X_temp;
                        end
                    end
                else
                    if obj.model_subtype==1
                        %X_z has more than one column
                        X_temp=[];
                        for i=1:length(obj.stem_varset_p.X_z)
                            X_temp2=[];
                            for j=1:size(obj.stem_varset_p.X_z{i},2)
                                X_temp2=cat(2,X_temp2,diag(obj.stem_varset_p.X_z{i}(:,j)));
                            end
                            X_temp=cat(1,X_temp,X_temp2);
                        end
                        X_z=X_temp;
                    else
                        %X_z has one column
                        X_z=[];
                        for i=1:length(obj.stem_varset_p.X_z)
                            X_z=cat(1,X_z,obj.stem_varset_p.X_z{i});
                        end
                    end
                end
                
                obj.X_z=X_z;
                clear X_z
            end

            %X_p
            if not(isempty(obj.stem_varset_p.X_p))
                X_p=[];
                for i=1:length(obj.stem_varset_p.X_p)
                    X_p=cat(1,X_p,obj.stem_varset_p.X_p{i});
                end
                obj.X_p=X_p;
                clear X_p;
            end
            disp('Generation ended.');
        end
        
        function update_M(obj)
            %DESCRIPTION: generates the vector M
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT         
            %
            %none: the vector M is generated ad updated
            
            disp('Generating M replication vector...');
            M=[];
            blocks=[0 cumsum(obj.stem_varset_b.dim)];
            for j=1:obj.stem_varset_p.nvar
                dmax=distdim(distance(0,0,obj.stem_gridlist_b.grid{j}.pixel_side_w,obj.stem_gridlist_b.grid{j}.pixel_side_w), obj.stem_gridlist_p.grid{1}.unit, 'km');
                for i=1:size(obj.stem_gridlist_p.grid{j}.coordinate,1)
                    d=distdim(distance(obj.stem_gridlist_p.grid{j}.coordinate(i,:),obj.stem_gridlist_b.grid{j}.coordinate), obj.stem_gridlist_p.grid{1}.unit, 'km');
                    [~,idx]=min(d);
                    if d>dmax
                        %warning(['Point ',num2str(i),' of point variable ',num2str(j),' does not belong to any pixel. The nearest pixel at ',num2str(m),' km is considered']);
                    end
                    M=cat(1,M,idx+blocks(j));
                end
            end
            obj.M=M;
            disp('Generation ended.');
        end
        
        function update_distance(obj,type,force)
            %DESCRIPTION: generates the distance matrices
            %
            %INPUT
            %obj    - [stem_data object] (1x1) the stem_data object
            %<type> - [string]           (1x1) (Default: 'both') 'point': only the distance matrix for the point data is evaluated. 
            %                                                    'pixel': only the distance matrix for the pixel data is evaluated.
            %                                                    'both':  both the matrices are evaluated.
            %
            %OUTPUT         
            %
            %none: the DistMat_p and DistMat_b property are generated and updated

            if nargin<2
                type='both';
            end
            if nargin<3
                force=0;
            end
            cmp=strcmp(type,{'both','point','pixel'});
            if sum(cmp)==0
                error('type must be point, pixel or both');
            end
            
            if strcmp(type,'point')||strcmp(type,'both')
                if not(isempty(obj.stem_varset_p.X_p))||force||obj.model_type==1
                    disp('Generating point distance matrices...');
                    obj.DistMat_p=obj.stem_gridlist_p.get_distance_matrix();
                    if obj.model_type==1
                        if obj.model_subtype==1
                            %y_t and the latent variable z_t have different dimensions
                            type=1;
                            idx_var=1;
                            obj.DistMat_z=kron(ones(size(obj.stem_varset_p.X_z,2))*obj.stem_gridlist_p.get_distance_matrix(type,idx_var));
                        end
                    end
                    disp('Generation ended.');
                end
            end
            if strcmp(type,'pixel')||strcmp(type,'both')
                if not(isempty(obj.stem_gridlist_b))&&not(isempty(obj.stem_varset_b.X_bp))
                    disp('Generating pixel data distance matrices...');
                    obj.DistMat_b=obj.stem_gridlist_b.get_distance_matrix(obj.pixel_correlated);
                    disp('Generation ended.');
                end
            end
        end     
       
        %Data transform
        function detrend_Y(obj)
            %DESCRIPTION: remove the mean from each time series in Y
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the Y property is updated
            
            disp('Point level data detrend started...');
            obj.stem_varset_p.detrend;
            disp('Point level data detrend ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data detrend started...');
                obj.stem_varset_b.detrend;
                disp('Pixel data detrend ended.');
            end
            disp('Updtaing data matrices after detrend...');
            obj.update_data;            
        end
        
        function standardize_Y(obj)
            %DESCRIPTION: each time series in Y is standardized
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the Y property is updated            
            
            disp('Point level data site by site standardization started...');
            obj.stem_varset_p.standardize_Y;
            disp('Point level data site by site standardization ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data site by site standardization started...');
                obj.stem_varset_b.standardize_Y;
                disp('Pixel data site by site standardization ended.');
            end
            disp('Updtaing data matrices after site by site standardization...');
            obj.update_data;
        end
        
        function standardize(obj)
            %DESCRIPTION: standardize the matrices Y, X_bp, X_beta, X_z and X_p with respect to their overall mean and overall standard deviation
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            disp('Point level data standardization started...');
            obj.stem_varset_p.standardize;
            disp('Point level data standardization ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data standardization started...');
                obj.stem_varset_b.standardize;
                disp('Pixel data standardization ended.');
            end
            disp('Updtaing data matrices after standardization...');
            obj.update_data;
        end
        
        function log_transform(obj)
            %DESCRIPTION: log-transforms the matrix Y
            %
            %INPUT
            %obj - [stem_data object] (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrix Y is updated            
            
            if obj.stem_varset_p.standardized
                error('The log_transform method must be called before the standardize method');
            end
            if not(isempty(obj.stem_varset_b))
                if obj.stem_varset_b.standardized
                    error('The log_transform method must be called before the standardize method');
                end
            end
            
            disp('Point level data log-transformation started...');
            obj.stem_varset_p.log_transform;
            disp('Point level data log-transformation ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data log-transformation started...');
                obj.stem_varset_b.log_transform;
                disp('Pixel data log-transformation ended.');
            end
            disp('Updtaing data matrices after log-transformation...');
            obj.update_data;
        end
        
        function time_average(obj,n_steps)
            %DESCRIPTION: computes time averages of n_steps for the matrice the matrices Y, X_bp, X_beta, X_z and X_p
            %
            %INPUT
            %obj        - [stem_data object] (1x1) the stem_data object
            %n_steps    - [integer >0]       (1x1) the number of temporal steps to average
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated
            
            if nargin<2
                error('The subsampling factor must be provided');
            end
            if n_steps<=1
                error('n_steps must be greater than 1');
            end
            if round(n_steps)~=n_steps
                error('n_steps must be an integer value');
            end
   
            obj.stem_varset_p.time_average(n_steps);
            if not(isempty(obj.stem_varset_b))
                obj.stem_varset_b.time_average(n_steps);   
            end
            indices=0:n_steps:obj.T;
            if indices(end)~=obj.T
                indices=[indices,obj.T];
            end
            obj.stem_datestamp.average_stamps(indices);
            
            obj.update_data;    
            disp('Time averaging ended.');
        end
        
        function time_crop(obj,dates_or_indices)
            %DESCRIPTION: crop the matrices Y, X_bp, X_beta, X_z and X_p with respect to time
            %
            %INPUT
            %obj                 - [stem_data object]       (1x1) the stem_data object
            %dates_or_indices    - [string | integer >0]    {2x1}|(dTx1) a {2x1} cell vector of starting and ending date in the format dd-mm-yyyy or a dTx1 vector of (possibly non consecutive) temporal indices
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated            

            if iscell(dates_or_indices)
                start_date=datenum(dates_or_indices{1},'dd-mm-yyyy');
                end_date=datenum(dates_or_indices{2},'dd-mm-yyyy');
                if start_date>end_date
                    error('Starting date must be lower or equal to ending date');
                end
                sd=find(obj.stem_datestamp.stamp==start_date);
                ed=find(obj.stem_datestamp.stamp==end_date);
                if isempty(sd)||isempty(ed)
                    error('Dates out of range');
                end
                indices=sd:ed;
            else
                if min(dates_or_indices)<0
                    error('Starting date must be higher than zero');
                end
                if max(dates_or_indices)>obj.T
                    error('Ending date must be lower or equal to the total number of temporal steps');
                end
                indices=dates_or_indices;
            end
            
            disp('Temporal cropping...');
            obj.stem_varset_p.time_crop(indices);
            if not(isempty(obj.stem_varset_b))
                obj.stem_varset_b.time_crop(indices);
            end
            
            obj.stem_datestamp.subset_stamps(indices);
            
            %looks for line of all missing for the sparse grids of the point data
            changed=0;
            for i=1:obj.stem_varset_p.nvar
                indices=sum(isnan(obj.stem_varset_p.Y{i}),2)==size(obj.stem_varset_p.Y{i},2);
                if sum(indices)>0
                    obj.site_crop('point',obj.stem_varset_p.Y_name{i},indices);
                    disp(['Deleted ',num2str(sum(indices)),' site(s) for the point variable ',obj.stem_varset_p.Y_name{i},' due to all missing.']);
                    changed=1;
                end
            end
            if changed
                disp('Updating point distance matrix after cropping...');
                obj.update_distance('point'); %only point because the pixel data are not deleted from the data matrix even if they are NaN for all th time steps
                disp('Update ended.');
                if not(isempty(obj.stem_varset_b))
                    disp('Updating M replication vector after cropping...');
                    obj.update_M;
                    disp('Update ended.');
                end
            end
            disp('Updating data matrix after cropping...');
            obj.update_data;
            disp('Update ended.');
            disp('Temporal cropping ended.');
        end     
        
        function space_crop(obj,box)
            %DESCRIPTION: crop the matrices Y, X_bp, X_beta, X_z and X_p with respect to space
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1) the stem_data object
            %box                 - [double]             (4x1) the geographic box [lat_min,lat_max,lon_min,lon_max]
            %
            %OUTPUT
            %
            %none: the matrices listed above are updated   
            
            lat_min=box(1);
            lat_max=box(2);
            lon_min=box(3);
            lon_max=box(4);
            if nargin<2
                error('The bounding-box must be provided');
            end
            if lat_min>lat_max
                error('The lat_min value must be lower than the lat_max value');
            end
            if lon_min>lon_max
                error('The lon_min value must be lower than the lon_max value');
            end
            disp('Point level data spatial cropping...');
            for i=1:obj.stem_varset_p.nvar
                GY=obj.stem_gridlist_p.grid{i}.coordinate;
                indices = GY(:,2)<lon_min | GY(:,2)>lon_max | GY(:,1)<lat_min | GY(:,1)>lat_max;
                if sum(indices)>0
                    if strcmp(obj.stem_gridlist_p.grid{i}.grid_type,'regular')
                        grid_lat=obj.stem_gridlist_p.grid{i}.coordinate(:,1);
                        grid_lon=obj.stem_gridlist_p.grid{i}.coordinate(:,2);
                        grid_lat=reshape(grid_lat,obj.stem_gridlist_p.grid{i}.grid_size);
                        grid_lat=grid_lat(:,1);
                        grid_lon=reshape(grid_lon,obj.stem_gridlist_p.grid{i}.grid_size);
                        grid_lon=grid_lon(1,:);
                        grid_lat(grid_lat<lat_min)=[];
                        grid_lat(grid_lat>lat_max)=[];
                        grid_lon(grid_lon<lon_min)=[];
                        grid_lon(grid_lon>lon_max)=[];
                    end
                    if strcmp(obj.stem_gridlist_p.grid{i}.grid_type,'regular')
                        obj.stem_gridlist_p.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                    end                    
                    obj.site_crop('point',obj.stem_varset_p.Y_name{i},indices);
                else
                    error(['    Variable ',obj.stem_varset_p.Y_name{i},' does not have sites in the crop area.']);
                end
            end
            disp('Point level data spatial cropping ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Pixel data spatial cropping...');
                for i=1:obj.stem_varset_b.nvar
                    GY=obj.stem_gridlist_b.grid{i}.coordinate;
                    indices = GY(:,2)<lon_min | GY(:,2)>lon_max | GY(:,1) < lat_min | GY(:,1) > lat_max;
                    if sum(indices)>0
                        if strcmp(obj.stem_gridlist_b.grid{i}.grid_type,'regular')
                            grid_lat=obj.stem_gridlist_b.grid{i}.coordinate(:,1);
                            grid_lon=obj.stem_gridlist_b.grid{i}.coordinate(:,2);
                            grid_lat=reshape(grid_lat,obj.stem_gridlist_b.grid{i}.grid_size);
                            grid_lat=grid_lat(:,1);
                            grid_lon=reshape(grid_lon,obj.stem_gridlist_b.grid{i}.grid_size);
                            grid_lon=grid_lon(1,:);
                            grid_lat(grid_lat<lat_min)=[];
                            grid_lat(grid_lat>lat_max)=[];
                            grid_lon(grid_lon<lon_min)=[];
                            grid_lon(grid_lon>lon_max)=[];
                        end
                        if strcmp(obj.stem_gridlist_b.grid{i}.grid_type,'regular')
                            obj.stem_gridlist_b.grid{i}.grid_size=[length(grid_lat),length(grid_lon)];
                        end
                        obj.site_crop('pixel',obj.stem_varset_b.Y_name{i},indices);
                    else
                        error(['    Variable ',obj.stem_varset_b.Y_name{i},' does not have sites in the crop area.']);
                    end
                end
                disp('Pixel data spatial cropping ended.');
            end
            disp('Updating data matrices after cropping...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after cropping...');
            obj.update_distance;
            disp('Update ended.');

            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after cropping...');
                obj.update_M;
                disp('Update ended.');
            end
        end   
        
        function site_crop(obj,type,var_name,indices,crossval)
            %DESCRIPTION: remove specific sites from the dataset
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1)  the stem_data object
            %type                - [string]             (1x1)  'point': remove the sites from the point dataset; 'pixel': remove the sites from the pixel dataset
            %var_name            - [string]             (1x1)  the name of the variable from which to remove the sites
            %indices             - [integer >0]         (dNx1) the indices of the sites to remove
            %crossval            - [integer]            (1x1)  1:sites are cropped for cross-validation, 0: sites are cropped to be removed
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated
            
            if sum(strcmp(type,{'point','pixel'}))==0
                error('Type must be either point or pixel');
            end
            if min(indices<1)
                error('The minimum value of indices cannot be lower than 1');
            end
            
            if nargin<5
                crossval=0;
            end
            
            if strcmp(type,'point')
                idx_var=obj.stem_varset_p.get_Y_index(var_name);
                
                if not(crossval)
                    %adjusting cross-validation sites before site cropping
                    if not(isempty(obj.stem_crossval))
                        for j=1:length(obj.stem_crossval.variable_name)
                            if (strcmp(obj.stem_crossval.variable_name{j},var_name))
                                disp(['Adjusting cross-validation site indices for point variable ',var_name,' due to site cropping.']);
                                L=false(obj.stem_varset_p.dim(idx_var),1);
                                L(obj.stem_crossval.indices{j})=true;
                                n_before=length(obj.stem_crossval.indices{j});
                                L(indices)=[];
                                obj.stem_crossval.indices{j}=find(L);
                                n_after=length(obj.stem_crossval.indices{j});
                                if n_after<n_before
                                    disp(['Removed ',num2str(n_before-n_after),' cross-validation sites for point variable ',var_name,' due to site cropping.']);
                                end
                            end
                        end
                    end
                end
                
                disp('Cropping sites...');
                obj.stem_varset_p.site_crop(var_name,indices);
                obj.stem_gridlist_p.grid{idx_var}.coordinate(indices,:)=[];
            else
                if isempty(obj.stem_varset_b)
                    error('No pixel data');
                end
                obj.stem_varset_b.site_crop(var_name,indices);
                idx_var=obj.stem_varset_b.get_Y_index(var_name);
                obj.stem_gridlist_b.grid{idx_var}.coordinate(indices,:)=[];
            end
            disp('Cropping ended.');
            disp('Updating data matrices after cropping...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after cropping...');
            obj.update_distance(type);
            disp('Update ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after cropping...');
                obj.update_M;
                disp('Update ended.');
            end
        end
        
        function missing_crop(obj,type,threshold)
            %DESCRIPTION: remove sites of the point variables ONLY with a missing data rate higher or equal to a threshold
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1) the stem_data object
            %type                - [string]             (1x1) 'point': remove the sites from the point dataset; 'pixel': remove the sites from the pixel dataset
            %threshold           - [double >0 and <1]   (1x1) the threshold
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated

            if sum(strcmp(type,{'point','pixel'}))==0
                error('Type must be either point or pixel');
            end
            if not(threshold>0 && threshold<1)
                error('The threshold must be >0 and <1');
            end
            
            disp('Cropping sites with missing data rate above threshold...');
            if strcmp(type,'point')
                idx_vec = obj.stem_varset_p.missing_crop(threshold);
                for i=1:length(idx_vec)
                   obj.stem_gridlist_p.grid{i}.coordinate(idx_vec{i},:)=[]; 
                end
            else
                if isempty(obj.stem_varset_b)
                    error('No pixel data');
                end
                idx_vec = obj.stem_varset_b.missing_crop(threshold);
                for i=1:length(idx_vec)
                    obj.stem_gridlist_b.grid{i}.coordinate(idx_vec{i},:)=[];
                end
            end
            
            disp('Cropping ended.');
            
            disp('Updating data matrices after cropping...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after cropping...');
            obj.update_distance(type);
            disp('Update ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after cropping...');
                obj.update_M;
                disp('Update ended.');
            end
        end
        
        function polygon_crop(obj,type,Xp,Yp,direction)
            %DESCRIPTION: remove sites inside or outside a polygon
            %
            %INPUT
            %obj                 - [stem_data object]   (1x1) the stem_data object
            %type                - [string]             (1x1) 'point': remove the sites from the point dataset; 'pixel': remove the sites from the pixel dataset
            %Xp                  - [double]             (Dx1) the X coordinates of the polygon (longitude if Xp is given in degrees)
            %Yp                  - [double]             (Dx1) the Y coordinates of the polyton (latitude if Yp is given in degrees)
            %direction           - [string]             (1x1) 'inside': remove the sites inside the polygon; 'outside': remove the sites outside the polygon
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated
            if nargin<5
                error('All the input arguments must be provided');
            end
            if sum(strcmp(type,{'point','pixel'}))==0
                error('Type must be either point or pixel');
            end
            if sum(strcmp(direction,{'inside','outside'}))==0
                error('direction must be either inside or outside');
            end
            if not(length(Xp)==length(Yp))
                error('Xp and Yp must have the same length');
            end
            if strcmp('outside',direction)
                disp('Cropping sites outside polygon...');
            else
                disp('Cropping sites inside polygon...');
            end
            if strcmp(type,'point')
                for i=1:obj.stem_varset_p.nvar
                    c=obj.stem_gridlist_p.grid{i}.coordinate;
                    IN=inpolygon(c(:,2),c(:,1),Xp,Yp);
                    if strcmp('outside',direction)
                        IN=not(IN);
                    end
                    idx=find(IN);
                    obj.site_crop('point',obj.stem_varset_p.Y_name{i},idx);
                end
            else
                if isempty(obj.stem_varset_b)
                    error('No pixel data');
                end
                for i=1:obj.stem_varset_b.nvar
                    c=obj.stem_gridlist_b.grid{i}.coordinate;
                    IN=inpolygon(c(:,2),c(:,1),Xp,Yp);
                    if strcmp('outside',direction)
                        IN=not(IN);
                    end
                    idx=find(IN);
                    obj.site_crop('pixel',obj.stem_varset_b.Y_name{i},idx);
                end
            end
            disp('Cropping ended.');
            
            disp('Updating data matrices after cropping...');
            obj.update_data;
            disp('Update ended.');
            disp('Updating distance matrices after cropping...');
            obj.update_distance(type);
            disp('Update ended.');
            if not(isempty(obj.stem_varset_b))
                disp('Updating M replication vector after cropping...');
                obj.update_M;
                disp('Update ended.');
            end
        end
        
        function remove_duplicated_sites(obj)
            %DESCRIPTION: remove the duplicated sites from each variable
            %
            %INPUT
            %obj        - [stem_data object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the matrices Y, X_bp, X_beta, X_z and X_p with are updated
            
            disp('Look for duplicated sites...')
            for i=1:obj.stem_varset_p.nvar
                idx=obj.stem_gridlist_p.grid{i}.duplicated_sites;
                if not(isempty(idx))
                    disp(['Removing ',num2str(length(idx)),' replicated sites for point variable ',obj.stem_varset_p.Y_name{i}]);
                    obj.site_crop('point',obj.stem_varset_p.Y_name{i},idx);
                    obj.stem_gridlist_p.grid{i}.duplicated_sites=[];
                end
            end
            
            if not(isempty(obj.stem_varset_b))
                for i=1:obj.stem_varset_b.nvar
                    idx=obj.stem_gridlist_b.grid{i}.duplicated_sites;
                    if not(isempty(idx))
                        disp(['Removing ',num2str(length(idx)),' replicated sites for pixel variable ',obj.stem_varset_b.Y_name{i}]);
                        obj.site_crop('pixel',obj.stem_varset_b.Y_name{i},idx);
                        obj.stem_gridlist_b.grid{i}.duplicated_sites=[];
                    end
                end
            end
            disp('Operation ended.')
        end        
        
        function h_fig = plot(obj,variable_name,variable_type,time_step,loading_name,loading_type)
            %DESCRIPTION: plot the variable data or the loading coefficients related to a variable
            %
            %INPUT
            %obj             - [stem_data object]   (1x1) the stem_data object
            %variable_name   - [string]             (1x1) the name of the variable
            %variable_type   - [string]             (1x1) the type of the variable, either 'point' or 'pixel'
            %time_step       - [integer >=0|string] (1x1) the time step to plot. If time_step=0 the temporal average is plotted. If time_step is a string it must be in the format dd-mm-yyyy
            %<loading_name>  - [string]             (1x1) (defalut: []) the name of the loading coefficients to plot (related to the variable specified)
            %<loading_type>  - [string]             (1x1) (default: []) the type of the loading coefficients. Can be 'beta', 'z', 'w_b' or 'w_p'
            %
            %OUTPUT
            %
            %h_fig           - [integer]            (1x1) the handle of the figure
            
            if nargin<4
                error('Not enough input arguments');
            end
            if nargin<5
                loading_name=[];
            end
            if nargin==5
                error('You must also provide the loading_type');
            end
            if ischar(time_step)
                date_num=datenum(time_step,'dd-mm-yyyy');
                time_step=find(obj.stem_datestamp.stamp==date_num);
                if isempty(time_step)
                    error(['The date stamp ',time_step,' cannot be found']);
                end
            end
            if time_step<0||time_step>obj.T
                error(['time_step out of bound. It must be between 0 and ',num2str(obj.T),' included']);
            end
            if not(strcmp(variable_type,'point')||strcmp(variable_type,'pixel'))
                error('variable_type can be either ''point'' or ''pixel''');
            end
            if nargin>=6
                if not(strcmp(loading_type,'X_beta')||strcmp(loading_type,'X_z')||strcmp(loading_type,'X_wp')||strcmp(loading_type,'X_wb'))
                    error('loading_type must be either ''X_beta'', ''X_z'', ''X_wb'' or ''X_wp''');
                end
                if strcmp(variable_type,'pixel')&&strcmp(loading_type,'X_wp')
                    error('The loading_type X_wp is not supported for pixel type data');
                end
            end
            
            if strcmp(variable_type,'point')
                index_var=obj.stem_varset_p.get_Y_index(variable_name);
                if isempty(index_var)
                    error(['The variable ',variable_name,' cannot be found as ',variable_type,' data']);
                end
                if not(isempty(loading_name))
                    if strcmp(loading_type,'X_beta')
                        indexl=obj.stem_varset_p.get_X_beta_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_beta for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_p.X_beta_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_beta{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_beta{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_z')
                        indexl=obj.stem_varset_p.get_X_z_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_z for the variable ',variable_name]);
                        end   
                        if not(obj.stem_varset_p.X_z_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_z{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_z{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_wb')
                        indexl=obj.stem_varset_p.get_X_bp_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_bp for the variable ',variable_name]);
                        end      
                        if not(obj.stem_varset_p.X_bp_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_bp{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_p.X_bp{index_var}(:,indexl,time_step);
                        end
                    end
                    if strcmp(loading_type,'X_wp')
                        indexl=obj.stem_varset_p.get_X_p_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_p for the variable ',variable_name]);
                        end     
                        if not(obj.stem_varset_p.X_p_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_p.X_p{index_var}(:,:,:,indexl),3));
                        else
                            data=obj.stem_varset_p.X_p{index_var}(:,:,time_step,indexl);
                        end
                    end
                else
                    if time_step==0
                        data=nanmean(obj.stem_varset_p.Y{index_var},2);
                    else
                        data=obj.stem_varset_p.Y{index_var}(:,time_step);
                    end
                end
                lat=obj.stem_gridlist_p.grid{index_var}.coordinate(:,1);
                lon=obj.stem_gridlist_p.grid{index_var}.coordinate(:,2);
                if isempty(loading_name)
                    if time_step==0
                        tit=['Average ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=[variable_name,' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                else
                    if time_step==0
                        tit=['Average loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=['Loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                end
                if strcmp(obj.stem_gridlist_p.grid{index_var}.unit,'deg')
                    xlab='Longitude';
                    ylab='Latitude';
                else
                    xlab=obj.stem_gridlist_p.grid{index_var}.unit;
                    ylab=obj.stem_gridlist_p.grid{index_var}.unit;
                end
                h=stem_misc.plot_map(lat,lon,data,obj.shape,tit,xlab,ylab);
            else
                index_var=obj.stem_varset_b.get_Y_index(variable_name);
                if isempty(index_var)
                    error(['The variable ',variable_name,' cannot be found as ',variable_type,' data']);
                end
                if not(isempty(loading_name))
                    if strcmp(loading_type,'X_beta')
                        indexl=obj.stem_varset_b.get_X_beta_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_beta for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_beta_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_beta{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_beta{index_var}(:,indexl,time_step);
                        end                        
                    end
                    if strcmp(loading_type,'X_z')
                        indexl=obj.stem_varset_b.get_X_z_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_z for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_z_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_z{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_z{index_var}(:,indexl,time_step);
                        end                        
                    end
                    if strcmp(loading_type,'X_wb')
                        indexl=obj.stem_varset_b.get_X_bp_index(loading_name,index_var);
                        if isempty(indexl)
                            error(['The loading coefficient ',loading_name,' cannot be found in X_bp for the variable ',variable_name]);
                        end
                        if not(obj.stem_varset_b.X_bp_tv)&&time_step>1
                            time_step=1;
                            disp('WARNING: the loading coefficient is time-invariant. t=1 is plotted');
                        end
                        if time_step==0
                            data=squeeze(nanmean(obj.stem_varset_b.X_bp{index_var}(:,indexl,:),3));
                        else
                            data=obj.stem_varset_b.X_bp{index_var}(:,indexl,time_step);
                        end                        
                    end
                else
                    if time_step==0
                        data=nanmean(obj.stem_varset_b.Y{index_var},2);
                    else
                        data=obj.stem_varset_b.Y{index_var}(:,time_step);
                    end
                end
                lat=obj.stem_gridlist_b.grid{index_var}.coordinate(:,1);
                lon=obj.stem_gridlist_b.grid{index_var}.coordinate(:,2);
                lat=reshape(lat,obj.stem_gridlist_b.grid{index_var}.grid_size);
                lon=reshape(lon,obj.stem_gridlist_b.grid{index_var}.grid_size);
                data=reshape(data,obj.stem_gridlist_b.grid{index_var}.grid_size);
                if isempty(loading_name)
                    if time_step==0
                        tit=['Average ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=[variable_name,' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                else
                    if time_step==0
                        tit=['Average loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' from ',datestr(obj.stem_datestamp.stamp(1)),' to ',datestr(obj.stem_datestamp.stamp(end))];
                    else
                        tit=['Loading coefficient ',loading_name,' in ',loading_type,' for variable ',variable_name, ' on ',datestr(obj.stem_datestamp.stamp(time_step))];
                    end
                end                
                if strcmp(obj.stem_gridlist_b.grid{index_var}.unit,'deg')
                    xlab='Longitude';
                    ylab='Latitude';
                else
                    xlab=obj.stem_gridlist_p.grid{index_var}.unit;
                    ylab=obj.stem_gridlist_p.grid{index_var}.unit;
                end                
                h=stem_misc.plot_map(lat,lon,data,obj.shape,tit,xlab,ylab);
            end
            if nargout>0
                h_fig=h;
            end
        end        
        
        function print(obj)
            %DESCRIPTION: print the information on data and their structure
            %
            %INPUT
            %obj  - [stem_data object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window
            
            stem_misc.disp_star('Data description');
            disp(['Number of variables: ',num2str(obj.nvar)]);
            disp(' ');
            disp('Point variables:');
            for i=1:length(obj.stem_varset_p.Y_name)
                disp(['  (',num2str(i),') - ',obj.stem_varset_p.Y_name{i}]);
            end
            disp(' '); 
            if not(isempty(obj.stem_varset_b))
                disp('Pixel variables:');
                for i=1:length(obj.stem_varset_b.Y_name)
                    disp(['  (',num2str(i),') - ',obj.stem_varset_b.Y_name{i}]);
                end
                disp(' ');
            end
            disp('Date and time steps:')
            disp(['  Stating date  : ',datestr(obj.stem_datestamp.date_start)]);
            disp(['  Ending date   : ',datestr(obj.stem_datestamp.date_end)]);
            disp(['  Temporal steps: ',num2str(obj.T)]);
            disp(' ');    
            
            disp('Bounding box of the point data:');
            if strcmp(obj.stem_gridlist_p.grid{1}.unit,'deg')
                prefix_x='  Longitude ';
                prefix_y='  Latitude ';
                postfix='°';
            else
                prefix_x='X ';
                prefix_y='Y ';
                postfix=[' ',obj.stem_gridlist_p.grid{1}.unit];
            end
            disp([prefix_y,'min : ',num2str(obj.stem_gridlist_p.box(1),'%05.2f'),postfix])
            disp([prefix_y,'max : ',num2str(obj.stem_gridlist_p.box(2),'%05.2f'),postfix])
            disp([prefix_x,'min: ',num2str(obj.stem_gridlist_p.box(3),'%05.2f'),postfix])
            disp([prefix_x,'max: ',num2str(obj.stem_gridlist_p.box(4),'%05.2f'),postfix])
            disp(' ');
            if not(isempty(obj.stem_varset_b))
                disp('Bounding box of the pixel data:');
                disp([prefix_y,'min : ',num2str(obj.stem_gridlist_b.box(1),'%05.2f'),postfix])
                disp([prefix_y,'max : ',num2str(obj.stem_gridlist_b.box(2),'%05.2f'),postfix])
                disp([prefix_x,'min: ',num2str(obj.stem_gridlist_b.box(3),'%05.2f'),postfix])
                disp([prefix_x,'max: ',num2str(obj.stem_gridlist_b.box(4),'%05.2f'),postfix])
                disp(' ');
            end

            disp('Variable description - Point data');
            output{1,1}='Name';
            output{1,2}='#sites';
            output{1,3}='Mean';            
            output{1,4}='Std';            
            output{1,5}='Min';            
            output{1,6}='Max';            
            output{1,7}='Missing';            
            for i=1:obj.stem_varset_p.nvar
                output{i+1,1}=obj.stem_varset_p.Y_name{i};
                output{i+1,2}=num2str(obj.stem_varset_p.dim(i),'%03.0f');
                output{i+1,3}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,4}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.Y{i})),'%05.2f');
                output{i+1,5}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,6}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.Y{i})),'%+05.2f');
                output{i+1,7}=[num2str(sum(isnan(stem_misc.vec(obj.stem_varset_p.Y{i})))/length(stem_misc.vec(obj.stem_varset_p.Y{i}))*100,'%05.2f'),'%'];
            end
            disp(output);
            output=[];
            if not(isempty(obj.stem_varset_b))
                disp('Variable description - Pixel data');
                output=cell(obj.stem_varset_b.nvar+1,7);
                output{1,1}='Name';
                output{1,2}='#sites';
                output{1,3}='Mean';
                output{1,4}='Std';
                output{1,5}='Min';
                output{1,6}='Max';
                output{1,7}='Missing';
                for i=1:obj.stem_varset_b.nvar
                    output{i+1,1}=obj.stem_varset_b.Y_name{i};
                    output{i+1,2}=num2str(obj.stem_varset_b.dim(i),'%03.0f');
                    output{i+1,3}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,4}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.Y{i})),'%05.2f');
                    output{i+1,5}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,6}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.Y{i})),'%+05.2f');
                    output{i+1,7}=[num2str(sum(isnan(stem_misc.vec(obj.stem_varset_b.Y{i})))/length(stem_misc.vec(obj.stem_varset_b.Y{i}))*100,'%05.2f'),'%'];
                end
            end
            disp(output);
            for i=1:obj.stem_varset_p.nvar
                stem_misc.disp_star(['Loading coefficients of the point variable ',obj.stem_varset_p.Y_name{i}]);
                if not(isempty(obj.stem_varset_p.X_beta_name))
                    disp('Loading coefficients related to beta:');
                    output=cell(length(obj.stem_varset_p.X_beta_name{i})+1,5);
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';
                    for j=1:length(obj.stem_varset_p.X_beta_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_beta_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_beta{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
                if not(isempty(obj.stem_varset_p.X_z_name))
                    disp('Loading coefficients related to the latent variable z(t):');
                    output=cell(length(obj.stem_varset_p.X_z_name{i})+1,5);
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                    
                    for j=1:length(obj.stem_varset_p.X_z_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_z_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_z{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
               if not(isempty(obj.stem_varset_p.X_p_name))
                    disp('Loading coefficients related to the latent variable w_p(s,t):');
                    output=cell(length(obj.stem_varset_p.X_p_name{i})+1,5);
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                       
                    for j=1:length(obj.stem_varset_p.X_p_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_p_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_p{i}(:,:,:,j))),'%+05.2f');
                    end
                    disp(output);
                end
                if not(isempty(obj.stem_varset_p.X_bp_name))
                    disp('Loading coefficients related to the latent variable w_b(s,t):');
                    output=cell(length(obj.stem_varset_p.X_bp_name{i})+1,5);
                    output{1,1}='Name';
                    output{1,2}='Mean';
                    output{1,3}='Std';
                    output{1,4}='Min';
                    output{1,5}='Max';                      
                    for j=1:length(obj.stem_varset_p.X_bp_name{i})
                        output{j+1,1}=obj.stem_varset_p.X_bp_name{i}{j};
                        output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                        output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%05.2f');
                        output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                        output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_p.X_bp{i}(:,j,:))),'%+05.2f');
                    end
                    disp(output);
                end
            end
                
            if not(isempty(obj.stem_varset_b))
                for i=1:obj.stem_varset_b.nvar
                    stem_misc.disp_star(['Loading coefficients of the pixel variable ',obj.stem_varset_b.Y_name{i}]);
                    if not(isempty(obj.stem_varset_b.X_beta_name))
                        disp('Loading coefficients related to beta:');
                        output=cell(length(obj.stem_varset_b.X_beta_name{i})+1,5);
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_beta_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_beta_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_beta{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                    if not(isempty(obj.stem_varset_b.X_z_name))
                        disp('Loading coefficients related to the latent variable z(t):');
                        output=cell(length(obj.stem_varset_b.X_z_name{i})+1,5);
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_z_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_z_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_z{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                    if not(isempty(obj.stem_varset_b.X_bp_name))
                        disp('Loading coefficients related to the latent variable w_b(s,t):');
                        output=cell(length(obj.stem_varset_b.X_bp_name{i})+1,5);
                        output{1,1}='Name';
                        output{1,2}='Mean';
                        output{1,3}='Std';
                        output{1,4}='Min';
                        output{1,5}='Max';
                        for j=1:length(obj.stem_varset_b.X_bp_name{i})
                            output{j+1,1}=obj.stem_varset_b.X_bp_name{i}{j};
                            output{j+1,2}=num2str(nanmean(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                            output{j+1,3}=num2str(nanstd(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%05.2f');
                            output{j+1,4}=num2str(nanmin(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                            output{j+1,5}=num2str(nanmax(stem_misc.vec(obj.stem_varset_b.X_bp{i}(:,j,:))),'%+05.2f');
                        end
                        disp(output);
                    end
                end
            end
        end
        
        %Export methods
        function N = N(obj)
            N=obj.stem_varset_p.N;
            if not(isempty(obj.stem_varset_b))
                N=N+obj.stem_varset_b.N;
            end
        end
        
        function Nb = Nb(obj)
            if not(isempty(obj.stem_varset_b))
                Nb = obj.stem_varset_b.N;
            else
                Nb=0;
            end
        end
        
        function Np = Np(obj)
            Np = obj.stem_varset_p.N;
        end
        
        function T = T(obj)
            T=obj.stem_varset_p.T;
        end
        
        function nvar = nvar(obj)
            nvar=obj.stem_varset_p.nvar;
            if not(isempty(obj.stem_varset_b))
                nvar=nvar+obj.stem_varset_b.nvar;
            end
        end
        
        function dim = dim(obj)
            dim=obj.stem_varset_p.dim;
            if not(isempty(obj.stem_varset_b))
                dim=[dim obj.stem_varset_b.dim];
            end
        end
        
        %Class set methods
        
        function standardized = get.standardized(obj)
            standardized=obj.stem_varset_p.standardized;
        end
        
        function log_transformed = get.log_transformed(obj)
            log_transformed=obj.stem_varset_p.log_transformed;
        end        
        
        function X_bp_tv = get.X_bp_tv(obj)
            if not(isempty(obj.X_bp))
                if size(obj.X_bp,3)==1
                    X_bp_tv=0;
                else
                    X_bp_tv=1;
                end
            else
                X_bp_tv=0;
            end
        end
        
        function X_beta_tv = get.X_beta_tv(obj)
            if not(isempty(obj.X_beta))
                if size(obj.X_beta,3)==1
                    X_beta_tv=0;
                else
                    X_beta_tv=1;
                end
            else
                X_beta_tv=0;
            end
        end
        
        function X_z_tv = get.X_z_tv(obj) 
            if not(isempty(obj.X_z))
                if size(obj.X_z,3)==1
                    X_z_tv=0;
                else
                    X_z_tv=1;
                end
            else
                X_z_tv=0;
            end
        end
        
        function X_p_tv = get.X_p_tv(obj)
            if not(isempty(obj.X_p))
                if size(obj.X_p,3)==1
                    X_p_tv=0;
                else
                    X_p_tv=1;
                end
            else
                X_p_tv=0;
            end
        end
   
        function X_tv = get.X_tv(obj)
            X_tv=obj.X_bp_tv | obj.X_z_tv | obj.X_p_tv;
        end
        
        function model_subtype = get.model_subtype(obj)
            if not(isempty(obj.stem_varset_p.X_z))
                if size(obj.stem_varset_p.X_z{1},2)==1
                    model_subtype=0;
                else
                    model_subtype=1;
                end
            else
                model_subtype=[];
            end
        end
        
        function set.stem_varset_p(obj,stem_varset_p)
           if not(isa(stem_varset_p,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           obj.stem_varset_p=stem_varset_p;
        end
        
        function set.stem_varset_b(obj,stem_varset_b)
           if not(isa(stem_varset_b,'stem_varset'))
               error('stem_varset must be of class stem_varset');
           end
           
           %            if not(isempty(stem_varset_b.X_beta))
           %                error('X_beta must be empty in stem_varset_b');
           %            end
           %
           %            if not(isempty(stem_varset_b.X_z))
           %                error('X_z must be empty in stem_varset_b');
           %            end

           if not(isempty(stem_varset_b.X_p))
               error('X_p must be empty in stem_varset_b');
           end
          
           obj.stem_varset_b=stem_varset_b;
        end        
        
        function set.stem_gridlist_p(obj,stem_gridlist_p)
            if not(isa(stem_gridlist_p,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            obj.stem_gridlist_p=stem_gridlist_p;
        end
        
        function set.stem_gridlist_b(obj,stem_gridlist_b)
            if not(isa(stem_gridlist_b,'stem_gridlist'))
                error('stem_gridlist must be of class stem_gridlist');
            end
            obj.stem_gridlist_b=stem_gridlist_b;
        end  

        function set.stem_crossval(obj,stem_crossval)
            if not(isa(stem_crossval,'stem_crossval'))
                error('stem_crossval must be of class stem_crossval');
            end
            obj.stem_crossval=stem_crossval;
        end
        
        function set.stem_datestamp(obj,stem_datestamp)
            if not(isa(stem_datestamp,'stem_datestamp'))
                error('stem_datestamp must be a stem_datestamp object');
            end
            obj.stem_datestamp=stem_datestamp;
        end
        
        function set.shape(obj,shape)
            obj.shape=shape;
        end
        
        function set.model_type(obj,model_type)
            if not(model_type>=0 && model_type<=3)
                error('model_type must be between 0 and 3 included');
            end
            obj.model_type=model_type;
        end
        
    end
end
        
        


