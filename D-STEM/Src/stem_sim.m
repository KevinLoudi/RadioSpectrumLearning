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

classdef stem_sim < handle
    
    properties
        stem_model=[];      %[stem_model object]    (1x1) stem_model object
        nan_rate=[0 0];     %[double >0 and <1]     (2x1) missing data rates for the point data and the pixel data
        nan_pattern_par=[]; %[double >0]            (2x1) missing data spatial-pattern parameters for point data and pixel data. The missing data spatial pattern is simulated by considering an exponential correlation function in the form exp(-d/theta) with d the euclidean distance between two sites.
    end
    
    methods
        function obj = stem_sim(stem_model)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_model  - [stem_model object] (1x1)
            %
            %OUTPUT
            %obj         - [stem_sim object]   (1x1)
            
            if nargin>=1
                if isa(stem_model,'stem_model')
                    obj.stem_model=stem_model;
                else
                    error('The input argument must be of class stem_model');
                end
            end
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            %DESCRIPTION: data simulation. Note that only the Y matrix is simulated. The loading vectors are taken from the stem_model object
            %
            %INPUT
            %obj                - [stem_sim object]      (1x1)
            %<nan_rate>         - [double >0 and <1]     (2x1) missing data rates for the point data and the pixel data
            %<nan_pattern_par>  - [double >0]            (2x1) missing data spatial-pattern parameters for point data and pixel data. The missing data spatial pattern is simulated by considering an exponential correlation function in the form exp(-d/theta) with d the euclidean distance between two sites.
            %
            %OUTPUT
            %none: the matrix Y of the stem_data object in the stem_model object is updated
            
            disp('Simulation started...');
            T=obj.stem_model.stem_data.T;
            N=obj.stem_model.stem_data.N;
            if nargin>=2
                obj.nan_rate=nan_rate;
            else
                obj.nan_rate=[];
            end
            if nargin>=3
                obj.nan_pattern_par=nan_pattern_par;
            else
                obj.nan_pattern_par=[];
            end
            
            obj.stem_model.stem_par=obj.stem_model.stem_par_initial;
            
            if not(isempty(obj.nan_rate))
                soglia_nan_p=norminv(1-obj.nan_rate(1)/2,0,1);
                nanmat_p=exp(-obj.stem_model.stem_data.DistMat_p./obj.nan_pattern_par(1));
            else
                nanmat_p=[];
            end
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))&&(not(isempty(obj.nan_rate)))
                soglia_nan_b=norminv(1-obj.nan_rate(2)/2,0,1);
                nanmat_b=exp(-obj.stem_model.stem_data.DistMat_b./obj.nan_pattern_par(2));
            else
                nanmat_b=[];
            end
            nancov_p=[];
            if not(isempty(nanmat_p))
                for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                    nancov_p=blkdiag(nancov_p,stem_misc.get_block(obj.stem_model.stem_data.stem_varset_p.dim,i,obj.stem_model.stem_data.stem_varset_p.dim,i,nanmat_p));
                end
            end
            nancov_b=[];
            if not(isempty(nanmat_b))
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    nancov_b=blkdiag(nancov_b,stem_misc.get_block(obj.stem_model.stem_data.stem_varset_b.dim,i,obj.stem_model.stem_data.stem_varset_b.dim,i,nanmat_b));
                end
            end            
            
            [sigma_eps,sigma_W_b,sigma_W_p,~,~,sigma_eta,G_tilde_diag,j_bp,j_p,j_z] = obj.stem_model.get_sigma();
            if obj.stem_model.stem_data.model_type==1
                s_eta=sigma_eta;
                r=length(G_tilde_diag);
                G=sparse(1:r,1:r,G_tilde_diag,r,r);
                mu0=zeros(r,1);
                sigma0=eye(r);
            else
                s_eta=obj.stem_model.stem_par.sigma_eta;
                G=obj.stem_model.stem_par.G;
                mu0=zeros(obj.stem_model.stem_par.p,1);
                sigma0=eye(obj.stem_model.stem_par.p);
            end
            
            
            if obj.stem_model.stem_par.p>0
                Z=stem_sim.ar1_sim(G,s_eta,T,mu0,sigma0);
            end
            
            Y=zeros(N,T);
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                if not(isempty(obj.stem_model.stem_data.stem_varset_b.X_bp))
                    W_b=mvnrnd(zeros(obj.stem_model.stem_data.stem_varset_b.N,1),sigma_W_b,T)';
                end
            end
            if obj.stem_model.stem_par.k>0
                W_p=zeros(obj.stem_model.stem_data.stem_varset_p.N,T,obj.stem_model.stem_par.k);
                for k=1:obj.stem_model.stem_par.k
                    W_p(:,:,k)=mvnrnd(zeros(obj.stem_model.stem_data.stem_varset_p.N,1),sigma_W_p{k},T)';
                end
            end
            W_eps=mvnrnd(zeros(N,1),sigma_eps,T)';
            
            for t=1:T
                if not(isempty(obj.stem_model.stem_data.X_bp))
                    if obj.stem_model.stem_data.X_bp_tv
                        Y(:,t)=Y(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(W_b(:,t),obj.stem_model.stem_data.M,'l'),obj.stem_model.stem_data.X_bp(:,1,t),'l'),j_bp,'l');
                    else
                        Y(:,t)=Y(:,t)+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(W_b(:,t),obj.stem_model.stem_data.M,'l'),obj.stem_model.stem_data.X_bp(:,1,1),'l'),j_bp,'l');
                    end
                end
                if not(isempty(obj.stem_model.stem_data.X_beta))
                    if obj.stem_model.stem_data.X_beta_tv
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_beta(:,:,t)*obj.stem_model.stem_par.beta;
                    else
                        Y(:,t)=Y(:,t)+obj.stem_model.stem_data.X_beta(:,:,1)*obj.stem_model.stem_par.beta;
                    end
                end
                if not(isempty(obj.stem_model.stem_data.X_z))
                    if obj.stem_model.stem_data.model_type==1
                        if obj.stem_model.stem_data.X_z_tv
                            Y(:,t)=Y(:,t)+diag(stem_misc.D_apply(obj.stem_model.stem_data.X_z(:,:,t),j_z,'l'))*Z(:,t);
                        else
                            Y(:,t)=Y(:,t)+diag(stem_misc.D_apply(obj.stem_model.stem_data.X_z(:,:,1),j_z,'l'))*Z(:,t);
                        end
                    else
                        if obj.stem_model.stem_data.X_z_tv
                            Y(:,t)=Y(:,t)+stem_misc.D_apply(obj.stem_model.stem_data.X_z(:,:,t),j_z,'l')*Z(:,t);
                        else
                            Y(:,t)=Y(:,t)+stem_misc.D_apply(obj.stem_model.stem_data.X_z(:,:,1),j_z,'l')*Z(:,t);
                        end
                    end
                end      
                if not(isempty(obj.stem_model.stem_data.X_p))
                    if obj.stem_model.stem_data.X_p_tv
                        for k=1:obj.stem_model.stem_par.k
                            Y(:,t)=Y(:,t)+stem_misc.D_apply(stem_misc.D_apply(W_p(:,t,k),obj.stem_model.stem_data.X_p(:,1,t,k),'l'),j_p(:,k),'l');
                        end
                    else
                        for k=1:obj.stem_model.stem_par.k
                            Y(:,t)=Y(:,t)+stem_misc.D_apply(stem_misc.D_apply(W_p(:,t,k),obj.stem_model.stem_data.X_p(:,1,1,k),'l'),j_p(:,k),'l');    
                        end
                    end
                end
                Y(:,t)=Y(:,t)+W_eps(:,t);

                if not(isempty(nancov_p)&&isempty(nancov_b))
                    nanfill_p=mvnrnd(zeros(size(nancov_p,1),1),nancov_p);
                    nanfill_p(abs(nanfill_p)>=soglia_nan_p)=NaN;
                    nanfill_p(abs(nanfill_p)<soglia_nan_p)=1;
                    if not(isempty(nancov_b))
                        nanfill_b=mvnrnd(zeros(size(nancov_b,1),1),nancov_b);
                        nanfill_b(abs(nanfill_b)>=soglia_nan_b)=NaN;
                        nanfill_b(abs(nanfill_b)<soglia_nan_b)=1;
                    else
                        nanfill_b=[];
                    end
                    nanfill=[nanfill_p,nanfill_b];
                    Y(:,t)=Y(:,t).*nanfill';
                end
            end

            blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_p.dim)];
            Y_temp=cell(obj.stem_model.stem_data.stem_varset_p.nvar,1);
            for i=1:obj.stem_model.stem_data.stem_varset_p.nvar
                Y_temp{i}=Y(blocks(i)+1:blocks(i+1),:);
            end
            obj.stem_model.stem_data.stem_varset_p.Y=Y_temp;
            if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                temp=max(blocks);
                blocks=[0 cumsum(obj.stem_model.stem_data.stem_varset_b.dim)]+temp;
                Y_temp=cell(obj.stem_model.stem_data.stem_varset_b.nvar,1);
                for i=1:obj.stem_model.stem_data.stem_varset_b.nvar
                    Y_temp{i}=Y(blocks(i)+1:blocks(i+1),:);
                end
                obj.stem_model.stem_data.stem_varset_b.Y=Y_temp;                
            end
            obj.stem_model.stem_data.update_data();
            obj.stem_model.stem_data.simulated=1;
            disp('Simulation ended. New Y updated.');
            disp('');
        end
        
        %Class set methods
        function set.nan_rate(obj,nan_rate)
            if not(isempty(nan_rate))
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    if not(length(nan_rate)==2)
                        error('nan_rate must be a 2x1 vector');
                    end
                    for i=1:2
                        if (nan_rate(i)<0)||(nan_rate(i)>=1)
                            error('The nan_rate elements must be in the interval [0,1)');
                        end
                    end
                else
                    if not(isscalar(nan_rate))
                        error('nan_rate must be a scalar');
                    end
                    if (nan_rate<0)||(nan_rate>=1)
                        error('The nan_rate must be in the interval [0,1)');
                    end
                end
            end
            obj.nan_rate=nan_rate;
        end
        
        function set.nan_pattern_par(obj,nan_pattern_par)
            if not(isempty(nan_pattern_par))
                if not(isempty(obj.stem_model.stem_data.stem_varset_b))
                    if not(length(nan_pattern_par)==2)
                        error('nan_pattern_par must be a 2x1 vector');
                    end
                    for i=1:2
                        if nan_pattern_par(i)<=0
                            error('The nan_pattern_par elements must be > 0');
                        end
                    end
                else
                    if not(isscalar(nan_pattern_par))
                        error('nan_rate must be a scalar');
                    end
                    if nan_pattern_par<=0
                        error('The nan_pattern_par must be between > 0');
                    end
                end
            end
            obj.nan_pattern_par=nan_pattern_par;
        end
    end
    
    methods (Static)
        function Z = ar1_sim(G,sigma_eta,T,mu0,sigma_0)
            Z=zeros(length(G),T+1);
            Z0=mvnrnd(mu0,sigma_0)';
            Z(:,1)=Z0;
            for t=2:T+1
                Z(:,t)=G*Z(:,t-1)+mvnrnd(zeros(length(G),1),sigma_eta)';
            end
            Z=Z(:,2:T+1);
        end
    end
    
end

