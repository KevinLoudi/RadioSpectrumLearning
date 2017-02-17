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

classdef stem_kalman < handle
    
    %CONSTANTS
    %N   = n1_p+...+nq_p+n1_b+...+nq_b - total number of observation sites
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N_b = n1_b+...+nq_b+n1_b+...+nq_b - total number of covariates
    %T   - number of temporal steps
    %TT = T if the space-time varying coefficients are time-variant and TT=1 if they are time-invariant    
    
    properties
        stem_model=[];  %[stem_model object] (1x1) stem_model object
    end
    
    methods
        function obj = stem_kalman(stem_model)
            %DESCRIPTION: constructor of the class stem_kalman
            %
            %INPUT
            %
            %stem_model      - [stem_model object]    (1x1) stem_model object
            %
            %OUTPUT
            %obj             - [stem_kalman object]   (1x1) stem_kalman object            
            if isa(stem_model,'stem_model')
                obj.stem_model=stem_model;
            else
                error('The input argument must be of class stem_model');
            end
        end
        
        function [st_kalmanfilter_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_eta,G_tilde_diag,sigma_geo,aj_bp,aj_p,aj_z,M] = filter(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel)
            %DESCRIPTION: Kalman filter front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]              (1x1)  stem_kalman object
            %<compute_logL>                 - [boolean]                         (1x1)  (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]                         (1x1)  (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]                      (dTx1) (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                          (1x1)  (defalut: []) full or relative path of the folder to use for distributed computation
            %    
            %OUTPUT
            %st_kalmanfilter_result         - [stem_kalmanfilter_result object] (1x1)     
            %sigma_eps                      - [double]                          (NxN) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_b                      - [double]                          (N_bxN_b) sigma_W_b matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_p                      - [double]                          {k}(N_px_Ng) the sigma_W_p matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                          (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_eta                      - [double]                          (r x r) variance-covariance matrix of eta when model_type=1
            %G_tilde_diag                   - [double]                          (r x 1) diagonal of the G_tilde matrix when model_type=1
            %sigma_geo                      - [double]                          (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_bp                          - [double]                          (Nx1) the aj_bp vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_p                           - [double]                          (Nx1) the aj_p vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_z                           - [double]                          (Nx1) the aj_z vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)     
            %M                              - [integer >0]                      (N_px1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
                
            if nargin<2
                compute_logL=0;
            end
            if nargin<3
                enable_varcov_computation=0;
            end
            if nargin<4
                pathparallel=[];
                time_steps=[];
            end
            if nargin==4
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman filter started...');
            ct1=clock;
           
            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;            
            
            [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
            
            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            tapering=obj.stem_model.tapering;  
            
            if obj.stem_model.stem_data.model_type==1
                s_eta=sigma_eta;
                r=length(G_tilde_diag);
                G=sparse(1:r,1:r,G_tilde_diag,r,r);
                z0=zeros(r,1);
                P0=eye(r);
            else
                s_eta=par.sigma_eta;
                G=par.G;
                z0=zeros(obj.stem_model.stem_par.p,1);
                P0=eye(obj.stem_model.stem_par.p);
            end
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter(data.Y,data.X_bp,data.X_beta,data.X_z,data.X_p,par.beta,G,s_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,obj.stem_model.stem_data.model_type,obj.stem_model.stem_data.model_subtype);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = stem_kalman.Kfilter_parallel(data.Y,data.X_bp,data.X_beta,data.X_z,data.X_p,par.beta,G,s_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,obj.stem_model.stem_data.model_type,obj.stem_model.stem_data.model_subtype);
            end
            st_kalmanfilter_result = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL);
            
            ct2=clock;
            disp(['    Kalman filter ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function [st_kalmansmoother_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_geo,aj_bp,aj_p,aj_z,M] = smoother(obj,compute_logL,enable_varcov_computation,time_steps,pathparallel)
            %DESCRIPTION: Kalman smoother front-end method
            %
            %INPUT
            %
            %obj                            - [stem_kalman object]    (1x1)  stem_kalman object
            %<compute_logL>                 - [boolean]               (1x1)  (default: 0) 1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %<enable_varcov_computation>    - [boolean]               (1x1)  (dafault: 0) 1:produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %<time_steps>                   - [integer >0]            (dTx1) (default: []) the subset of time steps with respect to which compute the Kalman filter
            %<pathparallel>                 - [string]                (1x1)  (defalut: []) full or relative path of the folder to use for distributed computation
            %    
            %OUTPUT
            %st_kalmansmoother_result       - [stem_kalmansmoother_result object] (1x1)     
            %sigma_eps                      - [double]                            (NxN) the sigma_eps matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_b                      - [double]                            (N_bxN_b) sigma_W_b matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_W_p                      - [double]                            {K}(N_pxN_p) the sigma_W_p matrices (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_Z                        - [double]                            (pxp) the sigma_Z matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %sigma_geo                      - [double]                            (NxN) the sigma_geo matrix (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_bp                          - [double]                            (Nx1) the aj_bp vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_p                           - [double]                            (Nx1) the aj_p vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            %aj_z                           - [double]                          (Nx1) the aj_z vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details) 
            %M                              - [integer >0]                        (N_px1) the M vector (passed as output to avoid recomputation. See the get_sigma method of the class stem_model for more details)
            
            if nargin<2
              compute_logL=0;
            end
            if nargin<3
                enable_varcov_computation=0;
            end
            if nargin<4
                pathparallel=[];
                time_steps=[];
            end
            if nargin==4
                error('The pathparallel input argument must be provided');
            end
            disp('    Kalman smoother started...');
            ct1=clock;

            data=obj.stem_model.stem_data;
            par=obj.stem_model.stem_par;
            
            [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = obj.stem_model.get_sigma();
            
            tapering=obj.stem_model.tapering;

            time_diagonal=obj.stem_model.stem_par.time_diagonal;
            
            if obj.stem_model.stem_data.model_type==1
                s_eta=sigma_eta;
                r=length(G_tilde_diag);
                G=sparse(1:r,1:r,G_tilde_diag,r,r);
                z0=zeros(r,1);
                P0=eye(r);
            else
                s_eta=par.sigma_eta;
                G=par.G;
                z0=zeros(obj.stem_model.stem_par.p,1);
                P0=eye(obj.stem_model.stem_par.p);
            end
            
            [zk_s,Pk_s,PPk_s,logL] = obj.Ksmoother(data.Y,data.X_bp,data.X_beta,data.X_z,...
                data.X_p,par.beta,G,s_eta,sigma_W_b,...
                sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,...
                time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,obj.stem_model.stem_data.model_type,obj.stem_model.stem_data.model_subtype);
            st_kalmansmoother_result = stem_kalmansmoother_result(zk_s,Pk_s,PPk_s,logL,obj.stem_model.stem_data.stem_datestamp);
            ct2=clock;
            disp(['    Kalman smoother ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
    end
    
    methods (Static)
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,model_type,model_subtype)
            %DESCRIPTION: Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     (Nx1xTT)    the full X_bp matrix
            %X_beta                         - [double]     (NxN_bxTT)  the full X_beta matrix
            %X_z                            - [double]     (NxpxTT)    the full X_z matrix
            %X_p                            - [double]     (Nx1xTTxK)  the full X_p matrix
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_type=1
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_type=1
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN)       variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %model_type                     - [integer >0] (1x1)       0: type 1 model, 1: type 2 model, 2: clustering model
            %model_subtype                  - [integer >0] (1x1)       currently used when model_type=1. 0: X_z{i} has only one column; 1: X_z{i} has more than one column
            % 
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)     the filtered state
            %zk_u                           - [double]     (pxT+1)     the updated state
            %Pk_f                           - [double]     (pxpxT+1)   variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     (pxpxT+1)   variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)       innovation vector at time t=T
            %J                              - [double]     (pxNxT+1)   innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)       observed-data log-likelihood
            
            if nargin<22
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
           
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1));
            if enable_varcov_computation
                J_all=zeros(p,size(Y,1),T+1);
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            
            logL=0;
            for t=2:T+1
                if size(X_z,3)==1
                    tK=2;
                else
                    tK=t; %time variant
                end
                
                if compute_sigma_geo
                    if not(isempty(X_bp))
                        sigma_geo=zeros(N);
                        if size(X_bp,3)>1
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,t-1),'b'),aj_bp,'b');
                        else
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,1),'b'),aj_bp,'b');
                        end
                    end

                    if not(isempty(X_p))
                        if isempty(X_bp)
                            if tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end                        
                        for k=1:size(X_p,4)
                            if size(X_p,3)>1
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,t-1,k),'b'),aj_p(:,k),'b');
                            else
                               sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,1,k),'b'),aj_p(:,k),'b');
                            end
                        end
                    end
                    if isempty(X_p)&&isempty(X_bp)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                end
                
                if size(X_beta,3)==1
                    tX=2;
                else
                    tX=t; %time variant
                end
                Lt=not(isnan(Y(:,t-1))); %note the t-1
                
                if (model_type==1)&&(model_subtype==0)
                    temp=X_z(:,:,tK-1);
                    temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                    X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                else
                    X_z_orlated=[X_z(:,:,tK-1);zeros(N-size(X_z(:,:,tK-1),1),size(X_z(:,:,tK-1),2))];
                end
                X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                
                X_z_orlated=X_z_orlated(Lt,:);
                if stem_misc.zero_density(X_z_orlated)>90
                    X_z_orlated=sparse(X_z_orlated);
                end
                
                X_beta_orlated=X_beta(:,:,tX-1);
                X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                X_beta_orlated=X_beta_orlated(Lt,:);
                if stem_misc.zero_density(X_beta_orlated)>90
                    X_beta_orlated=sparse(X_beta_orlated);
                end                
                
                temp=sigma_geo(Lt,Lt);
                if tapering
                    if not(stem_misc.isdiagonal(temp))
                        if compute_logL
                            r = symamd(temp);
                            c=chol(temp(r,r));
                            temp2=speye(sum(Lt));
                            temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
                            sigma_geo_inv=zeros(size(temp3));
                            sigma_geo_inv(r,:)=temp3;
                            clear temp2
                            clear temp3
                        end
                        temp=X_z_orlated'/temp;
                    else
                        d=1./diag(temp);
                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                        temp=X_z_orlated'*sigma_geo_inv;
                    end
                else
                    if not(stem_misc.isdiagonal(temp))
                        if compute_logL
                            c=chol(sigma_geo(Lt,Lt));
                            sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                        end
                        temp=X_z_orlated'/temp;
                    else
                        d=1./diag(temp);
                        sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                        temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                        for i=1:size(temp,1)
                            temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                        end
                    end
                end
                temp2=temp*X_z_orlated;
                
                if not(time_diagonal)
                    %FILTERING
                    zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                    
                    %UPDATING
                    
                    %Original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    %temp=X_z_orlated'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    %temp2=temp*X_z_orlated;
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                    end
                    
                    temp3=sparse(Pk_f(:,:,t)*X_z_orlated');
                    J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                         
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1); 
                    Pk_u(:,:,t)=(eye(p)-J(:,Lt)*X_z_orlated)*Pk_f(:,:,t);  %(6.22) Stoffer
                else
                    %FILTERING
                    zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                    Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                    
                    %UPDATING

                    %Original formula
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                    %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                    %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                    
                    %temp=X_z_orlated'*sigma_geo_inv; %note that temp can be computed in a distributed way before the KF is started
                    %temp2=temp*X_z_orlated;
                    P=diag(1./diag(Pk_f(:,:,t)));
                    if compute_logL
                        sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                    end
                    temp3=Pk_f(:,:,t)*X_z_orlated';
                    J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(P+temp2)*temp;
                    
                    if not(isempty(X_beta))
                        innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    else
                        innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                    end
                    
                    zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                    Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                end
                if compute_logL
                    if tapering
                        r = symamd(sigma_t_inv);
                        c=chol(sigma_t_inv(r,r));
                    else
                        c=chol(sigma_t_inv);
                    end
                    logL=logL+(-(2*sum(log(diag(c))))); %the negative sign is due to the fact that is the log of sigma_t that must be computed
                    logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                end
                clear temp
                clear temp2
                clear temp3
                if enable_varcov_computation
                    J_all(:,:,t)=J;
                end
            end
            logL=-logL/2;
            J_last=J;
            if enable_varcov_computation
                J=J_all;
            else
                J=[];
            end
        end
        
        function [zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL] = Kfilter_parallel(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,model_type,model_subtype)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     (Nx1xTT)    the full X_bp matrix
            %X_beta                         - [double]     (NxN_bxTT)  the full X_beta matrix
            %X_z                            - [double]     (NxpxTT)    the full X_z matrix
            %X_p                            - [double]     (Nx1xTTxK)  the full X_p matrix
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_type=1
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_type=1
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN)       variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)      time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)       full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %model_type                     - [integer >0] (1x1)       0: type 1 model, 1: type 2 model, 2: clustering model
            %model_subtype                  - [integer >0] (1x1)       currently used when model_type=1. 0: X_z{i} has only one column; 1: X_z{i} has more than one column
            %  
            %OUTPUT 
            %zk_f                           - [double]     (pxT+1)     the filtered state
            %zk_u                           - [double]     (pxT+1)     the updated state
            %Pk_f                           - [double]     (pxpxT+1)   variance-covariance matrix of the filtered state
            %Pk_u                           - [double]     (pxpxT+1)   variance-covariance matrix of the updated state
            %J_last                         - [double]     (pxN)       innovation vector at time t=T
            %J                              - [double]     (pxNxT+1)   innovation vector from time t=0 to time t=T
            %logL                           - [double]     (1x1)       observed-data log-likelihood
            
            if nargin<22
                error('You have to provide all the input arguments');
            end
                        
            if size(X_beta,2)~=length(beta)
                error('X_beta and beta are not conformable');
            end
           
            if size(G,1)~=size(G,2)
                error('G must be square');
            end
            
            if size(sigma_eta,1)~=size(sigma_eta,2)
                error('sigma_eta must be square');
            end      
            
            if size(G,1)~=size(sigma_eta,1)
                error('G and sigma_eta must have the same dimensions');
            end
            
            if size(Y,1)~=size(sigma_eps,1)
                error('The dimensions of sigma_eps must be equal to the number of rows of Y');
            end
            
            if size(z0,1)~=size(G,1)
                error('The length of z0 must be equal to the dimensions of G');
            end
            
            if size(P0,1)~=size(P0,2)
                error('P0 must be square');
            end
            
            if size(P0,1)~=size(G,1)
                error('The dimensions of P0 must be equal to the dimensions of G and sigma_eps');
            end
            
            if isempty(sigma_geo)
                compute_sigma_geo=1;
            else
                compute_sigma_geo=0;
            end
            
            time_steps=time_steps+1; %!!!
            
            min_ts=min(time_steps);
            max_ts=max(time_steps);
            if min_ts==2 %note the 2 due to +1
                server=1;
                %if min_ts==2 it means the call of the function is on the server
            else
                server=0;
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            zk_f=zeros(p,T+1);
            zk_u=zeros(p,T+1);
            Pk_f=zeros(p,p,T+1);
            Pk_u=zeros(p,p,T+1);
            J=zeros(p,size(Y,1));
            if enable_varcov_computation
                J_all=zeros(p,size(Y,1),T+1);
            end
            innovation=zeros(size(Y,1),1);
            
            zk_u(:,1)=z0;
            Pk_u(:,:,1)=P0;
            logL=0;
         
            if server
                for t=2:T+1
                    if size(X_z,3)==1
                        tK=2;
                    else
                        tK=t; 
                    end
                    
                    if compute_sigma_geo
                        if not(isempty(X_bp))
                            sigma_geo=zeros(N);
                            if size(X_bp,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,t-1),'b'),aj_bp,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,1),'b'),aj_bp,'b');
                            end
                        end
                        
                        if not(isempty(X_p))
                            if isempty(X_bp)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_p,4)
                                if size(X_p,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,t-1,k),'b'),aj_p(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,1,k),'b'),aj_p(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_p)&&isempty(X_bp)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end
                    
                    if size(X_beta,3)==1
                        tX=2;
                    else
                        tX=t; %time variant
                    end
                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                           
                    if (model_type==1)&&(model_subtype==0)
                        temp=X_z(:,:,tK-1);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[X_z(:,:,tK-1);zeros(N-size(X_z(:,:,tK-1),1),size(X_z(:,:,tK-1),2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    
                    X_z_orlated=X_z_orlated(Lt,:);
                    if stem_misc.zero_density(X_z_orlated)>90
                        X_z_orlated=sparse(X_z_orlated);
                    end
                    
                    X_beta_orlated=X_beta(:,:,tX-1);
                    X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                    X_beta_orlated=X_beta_orlated(Lt,:);
                    if stem_misc.zero_density(X_beta_orlated)>90
                        X_beta_orlated=sparse(X_beta_orlated);
                    end
                    
                    if t>max_ts
                        %wait for the proper file from the clients
                        exit=0;
                        %disp(['        Waiting for kalman_output_',num2str(t)]);
                        while not(exit)
                            exit=exist([pathparallel,'kalman_ouput_',num2str(t),'.mat'],'file');
                        end
                        read=0;
                        while not(read)
                            try
                                load([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                read=1;
                                %disp(['        kalman_ouput_',num2str(t),' readed']);
                            catch
                            end
                            pause(0.05);
                        end
                        deleted=0;
                        while not(deleted)
                            try
                                delete([pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                                deleted=1;
                                %disp(['        kalman_ouput_',num2str(t),' deleted']);
                            catch
                            end
                            pause(0.05);
                        end
                    end
                    
                    if not(time_diagonal)
                        %FILTERING
                        zk_f(:,t)=G*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=G*Pk_u(:,:,t-1)*G'+sigma_eta; %(6.20) Stoffer
                        
                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        %UPDATING
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=sigma_geo(Lt,Lt);
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
                                        sigma_geo_inv=zeros(size(temp3));
                                        sigma_geo_inv(r,:)=temp3;
                                        clear temp2
                                        clear temp3
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    temp=X_z_orlated'*sigma_geo_inv;
                                end
                            else
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                                    for i=1:size(temp,1)
                                        temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                                    end
                                end
                            end
                            temp2=temp*X_z_orlated;
                            temp3=Pk_f(:,:,t)*X_z_orlated';
                            J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                        else
                            %temp and temp2 has already been read from the file
                            temp3=Pk_f(:,:,t)*X_z_orlated';
                            J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(Pk_f(:,:,t)\eye(size(temp2))+temp2)*temp;
                            if compute_logL
                                temp_s=sigma_geo(Lt,Lt);
                                if tapering
                                    if not(stem_misc.isdiagonal(temp))
                                        r = symamd(temp_s);
                                        c=chol(temp_s(r,r));
                                        temp2_s=speye(sum(Lt));
                                        temp3_s=full(stem_misc.chol_solve(full(c),temp2_s(r,:)));
                                        sigma_geo_inv=zeros(size(temp3_s));
                                        sigma_geo_inv(r,:)=temp3_s;
                                        clear temp2_s
                                        clear temp3_s
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                else
                                    if not(stem_misc.isdiagonal(temp))
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                end
                            end
                        end
                        
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/((Pk_f(:,:,t)\eye(size(temp2)))+temp2))*temp;
                            if tapering
                                r = symamd(sigma_t_inv);
                                c=chol(sigma_t_inv(r,r));
                            else
                                c=chol(sigma_t_inv);
                            end
                            logL=logL+1/(2*sum(log(diag(c))));
                            logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                        end

                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u(:,:,t)=(eye(p)-J(:,Lt)*X_z_orlated)*Pk_f(:,:,t);  %(6.22) Stoffer
                    else
                        %FILTERING
                        zk_f(:,t)=diag(G).*zk_u(:,t-1); %(6.19) Stoffer
                        Pk_f(:,:,t)=diag(diag(G).^2.*diag(Pk_u(:,:,t-1))+diag(sigma_eta)); %(6.20) Stoffer
                        
                        %UPDATING
                        
                        %Original formula
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'/(X_z(Lt,i,tK-1)*Pk_f(i,i,t)*X_z(Lt,i,tK-1)'+sigma_geo(Lt,Lt)); %(6.23) Stoffer
                        %Sherman-Morrison-Woodbury formula: (B*P*B+D)^-1=D^-1-D^-1*B(P^-1+B*D^-1*B)^-1*B*D^-1
                        %J(i,Lt,t)=Pk_f(i,i,t)*X_z(Lt,i,tK-1)'*(sigma_geo_inv-sigma_geo_inv*X_z(Lt,i,tK-1)/(1/Pk_f(i,i,t)+X_z(Lt,i,tK-1)'*sigma_geo_inv*X_z(Lt,i,tK-1))*(X_z(Lt,i,tK-1)'*sigma_geo_inv));
                        
                        if not(isempty(X_beta))
                            innovation(Lt,1)=Y(Lt,t-1)-X_beta_orlated*beta-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        else
                            innovation(Lt,1)=Y(Lt,t-1)-X_z_orlated*zk_f(:,t); %(6.21) Stoffer %note the t-1 on Y and X
                        end
                        
                        if t<=max_ts %the time steps up to max_ts are computed locally
                            temp=sigma_geo(Lt,Lt);
                            if tapering
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        r = symamd(temp);
                                        c=chol(temp(r,r));
                                        temp2=speye(sum(Lt));
                                        temp3=full(stem_misc.chol_solve(full(c),temp2(r,:)));
                                        sigma_geo_inv=zeros(size(temp3));
                                        sigma_geo_inv(r,:)=temp3;
                                        clear temp2
                                        clear temp3
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    temp=X_z_orlated'*sigma_geo_inv;
                                end
                            else
                                if not(stem_misc.isdiagonal(temp))
                                    if compute_logL
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    end
                                    temp=X_z_orlated'/temp;
                                else
                                    d=1./diag(temp);
                                    sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                                    temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                                    for i=1:size(temp,1)
                                        temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                                    end
                                end
                            end
                            temp2=temp*X_z_orlated;
                        else
                            %temp and temp2 has already been read from the file
                            if compute_logL
                                temp_s=sigma_geo(Lt,Lt);
                                if tapering
                                    if not(stem_misc.isdiagonal(temp))
                                        r = symamd(temp_s);
                                        c=chol(temp_s(r,r));
                                        temp2_s=speye(sum(Lt));
                                        temp3_s=full(stem_misc.chol_solve(full(c),temp2_s(r,:)));
                                        sigma_geo_inv=zeros(size(temp3_s));
                                        sigma_geo_inv(r,:)=temp3_s;
                                        clear temp2_s
                                        clear temp3_s
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                else
                                    if not(stem_misc.isdiagonal(temp))
                                        c=chol(sigma_geo(Lt,Lt));
                                        sigma_geo_inv=stem_misc.chol_solve(c,eye(sum(Lt)));
                                    else
                                        d=1./diag(temp_s);
                                        sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                                    end
                                end
                            end
                        end
                        
                        P=diag(1./diag(Pk_f(:,:,t)));
                        if compute_logL
                            sigma_t_inv=sigma_geo_inv-(temp'/(P+temp2))*temp;
                            if tapering
                                r = symamd(sigma_t_inv);
                                c=chol(sigma_t_inv(r,r));
                            else
                                c=chol(sigma_t_inv);
                            end
                            logL=logL+1/(2*sum(log(diag(c))));
                            logL=logL+innovation(Lt,1)'*sigma_t_inv*innovation(Lt,1);
                        end
                        
                        temp3=Pk_f(:,:,t)*X_z_orlated';
                        J(:,Lt)=Pk_f(:,:,t)*temp-temp3*temp'/(P+temp2)*temp;      
                        
                        zk_u(:,t)=zk_f(:,t)+J(:,Lt)*innovation(Lt,1);
                        Pk_u(:,:,t)=diag(diag((eye(p)-J(:,Lt)*X_z_orlated)).*diag(Pk_f(:,:,t))); %(6.22) Stoffer
                    end
                    if enable_varcov_computation
                        J_all(:,:,t)=J;
                    end
                end
                logL=-logL/2;
                J_last=J;
                if enable_varcov_computation
                    J=J_all;
                else
                    J=[];
                end
            else
                %client computation
                for t=time_steps
                    if size(X_z,3)==1
                        tK=2;
                    else
                        tK=t;
                    end
                    if compute_sigma_geo
                        if not(isempty(X_bp))
                            sigma_geo=zeros(N);
                            if size(X_bp,3)>1
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,t-1),'b'),aj_bp,'b');
                            else
                                sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),X_bp(:,1,1),'b'),aj_bp,'b');
                            end
                        end
                        
                        if not(isempty(X_p))
                            if isempty(X_bp)
                                if tapering
                                    sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                                else
                                    sigma_geo=zeros(N);
                                end
                            end
                            for k=1:size(X_p,4)
                                if size(X_p,3)>1
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,t-1,k),'b'),aj_p(:,k),'b');
                                else
                                    sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},X_p(:,1,1,k),'b'),aj_p(:,k),'b');
                                end
                            end
                        end
                        if isempty(X_p)&&isempty(X_bp)
                            sigma_geo=sigma_eps;
                        else
                            sigma_geo=sigma_geo+sigma_eps;
                        end
                    end

                    Lt=not(isnan(Y(:,t-1))); %note the t-1
                    
                    if (model_type==1)&&(model_subtype==0)
                        temp=X_z(:,:,tK-1);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[X_z(:,:,tK-1);zeros(N-size(X_z(:,:,tK-1),1),size(X_z(:,:,tK-1),2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');

                    X_z_orlated=X_z_orlated(Lt,:);
                    if stem_misc.zero_density(X_z_orlated)>90
                        X_z_orlated=sparse(X_z_orlated);
                    end
                    
                    temp=sigma_geo(Lt,Lt);
                    if tapering
                        if not(stem_misc.isdiagonal(temp))
                            temp=X_z_orlated'/temp;
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d);
                            temp=X_z_orlated'*sigma_geo_inv;
                        end
                    else
                        if not(stem_misc.isdiagonal(temp))
                            temp=X_z_orlated'/temp;
                        else
                            d=1./diag(temp);
                            sigma_geo_inv=sparse(1:length(d),1:length(d),d); %sigma_geo_inv=diag(1./diag(temp));
                            temp=zeros(size(X_z_orlated,2),size(X_z_orlated,1));
                            for i=1:size(temp,1)
                                temp(i,:)=X_z_orlated(:,i).*diag(sigma_geo_inv);
                            end
                        end
                    end
                    temp2=temp*X_z_orlated;           
                    save([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],'temp','temp2');
                    movefile([pathparallel,'temp/kalman_ouput_',num2str(t),'.mat'],[pathparallel,'kalman_ouput_',num2str(t),'.mat']);
                end
                J_last=[];
                J=[];
            end
        end
        
        function [zk_s,Pk_s,PPk_s,logL] = Ksmoother(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,model_type,model_subtype)
            %DESCRIPTION: distributed Kalman filter implementation
            %
            %INPUT
            %
            %Y                              - [double]     (NxT)       the full observation matrix
            %X_bp                           - [double]     (Nx1xTT)    the full X_bp matrix
            %X_beta                         - [double]     (NxN_bxTT)  the full X_beta matrix
            %X_z                            - [double]     (NxpxTT)    the full X_z matrix
            %X_p                            - [double]     (Nx1xTTxK)  the full X_p matrix
            %beta                           - [double]     (N_bx1)     the beta model parameter
            %G                              - [double]     (pxp)|(rxr) the G model parameter or the G_tilde matrix when model_type=1
            %sigma_eta                      - [double]     (pxp)|(rxr) the sigma_eta model parameter or the sigma_eta matrix when model_type=1
            %sigma_W_b                      - [double]     (N_bxN_b)   variance-covariance matrix of W_b
            %sigma_W_p                      - [double]     {K}(N_pxN_p)variance-covariance matrices of the K W_p_i
            %sigma_eps                      - [double]     (NxN)       variance-covariance matrix of epsilon
            %sigma_geo                      - [double]     (NxN)       variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %aj_bp                          - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_p                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %aj_z                           - [double]     (Nx1)       see the details of the method get_aj of the class stem_model;
            %M                              - [integer >0] (N_px1)     see the details of the method update_M of the class stem_data            
            %z0                             - [double]     (px1)       the value of z at time t=0
            %P0                             - [double]     (pxp)       the variance-covariance matrix of z at time t=0
            %time_diagonal                  - [boolean]    (1x1)       1: G and sigma_eta are diagonal matrice; 0:otherwise
            %time_steps                     - [integer >0] (dTx1)      time steps with respect to which compute the Kalman filter
            %pathparallel                   - [string]     (1x1)       full or relative path of the folder to use for distributed computation
            %tapering                       - [boolean]    (1x1)       1: tapering is enabled; 0: tapering is not enabled
            %compute_logL                   - [boolean]    (1x1)       1: compute the observed-data log-likelihood; 0: the log-likelihood is not computed
            %enable_varcov_computation      - [boolean]    (1x1)       1: produce the output necessary to the computation of the variance-covariance matrix of the estimated model parameter; 0: the output is not produced
            %model_type                     - [integer >0] (1x1)       0: type 1 model, 1: type 2 model, 2: clustering model
            %model_subtype                  - [integer >0] (1x1)       currently used when model_type=1. 0: X_z{i} has only one column; 1: X_z{i} has more than one column
            % 
            %OUTPUT 
            %zk_s                           - [double]     (pxT+1)     the smoothed state
            %Pk_s                           - [double]     (pxpxT+1)   variance-covariance matrix of the smoothed state
            %PPk_s                          - [double]     (pxpxT+1)   lag-one variance-covariance matrix of the smoothed state
            %logL                           - [double]     (1x1)       observed-data log-likelihood
        
            if isempty(pathparallel)
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,tapering,compute_logL,enable_varcov_computation,model_type,model_subtype);
            else
                [zk_f,zk_u,Pk_f,Pk_u,J_last,~,logL] = stem_kalman.Kfilter_parallel(Y,X_bp,X_beta,X_z,X_p,beta,G,sigma_eta,sigma_W_b,sigma_W_p,sigma_eps,sigma_geo,aj_bp,aj_p,aj_z,M,z0,P0,time_diagonal,time_steps,pathparallel,tapering,compute_logL,enable_varcov_computation,model_type,model_subtype);
            end
            
            p=size(G,1);
            N=size(Y,1);
            T=size(Y,2);
            
            H=zeros(p,p,T+1);
            Pk_s=zeros(p,p,T+1);
            zk_s=zeros(p,T+1);
            zk_s(:,end)=zk_u(:,end); %inizializzazione (6.47) Stoffer
            Pk_s(:,:,end)=Pk_u(:,:,end); %inizializzazione (6.48) Stoffer
            PPk_s=zeros(p,p,T+1);
            
            for t=T+1:-1:2
                H(:,:,t-1)=Pk_u(:,:,t-1)*G'/(Pk_f(:,:,t)); %(6.49) Stoffer
                zk_s(:,t-1)=zk_u(:,t-1)+H(:,:,t-1)*(zk_s(:,t)-zk_f(:,t)); %(6.47) Stoffer
                Pk_s(:,:,t-1)=Pk_u(:,:,t-1)+H(:,:,t-1)*(Pk_s(:,:,t)-Pk_f(:,:,t))*H(:,:,t-1)'; %(6.48) Stoffer
            end
            
            Lt=not(isnan(Y(:,end)));
            
            if (model_type==1)&&(model_subtype==0)
                temp=X_z(:,:,end);
                temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
            else
                X_z_orlated=[X_z(:,:,end);zeros(N-size(X_z(:,:,end),1),size(X_z(:,:,end),2))];
            end
            X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');
            
            X_z_orlated=X_z_orlated(Lt,:);
            if stem_misc.zero_density(X_z_orlated)>90
                X_z_orlated=sparse(X_z_orlated);
            end
            
            PPk_s(:,:,end)=(eye(p)-J_last(:,Lt)*X_z_orlated)*G*Pk_u(:,:,end-1); %(6.55) Stoffer
            for t=T+1:-1:3
                PPk_s(:,:,t-1)=Pk_u(:,:,t-1)*H(:,:,t-2)'+H(:,:,t-1)*(PPk_s(:,:,t)-G*Pk_u(:,:,t-1))*H(:,:,t-2)'; %(6.56) Stoffer
            end
        end
    end
end