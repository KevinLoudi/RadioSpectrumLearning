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

classdef stem_model < handle
    
    %CONSTANTS
    %N   = n1_p+...+nq_p+n1_b+...+nq_b - total number of observation sites
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %K - the number of loading vectors related to the latent variable w_p
    %r - total number of elements of the latent variable z when model_type=1
    %S = 1 if only point data are considered and S=2 if both point and pixel data are considered
    
    properties
        stem_data=[];           %[stem_data object] (1x1) object containing all the data used to estimated the model
        stem_par=[];            %[stem_par object]  (1x1) parameter set updated at each iteration of the EM algorithm
        stem_par_initial=[];    %[stem_par object]  (1x1) starting parameter set
        stem_par_sim=[];        %[stem_par object]  (1x1) parameter set used to simulate data (if data are simulated)
        estimated=0;            %[boolean] (1x1) 0: the model has not been estimated; 1: the model has been estimated
    end
    
    properties (SetAccess = private)
        stem_EM_result=[];      %[stem_EM_result object] (1x1) object containing all the results of the EM estimation
        cross_validation=0;     %[boolean] (1x1) 0: the model has been estimated considering all the data; 1: the model has bee estimated excluding the cross-validation data.
        system_size=100;        %[integer] (1x1) if N > system_size than only the diagonal is computed in matrix multiply operations
        tapering=[];            %[boolean] (1x1) 0:tapering is not enabled; 1:tapering is enabled on point sites or pixel sites
        tapering_b=[];          %[boolean] (1x1) 0:tapering is not enabled on pixel sites; 1:tapering is enabled on pixel sites
        tapering_p=[];          %[boolean] (1x1) 0:tapering is not enabled on point sites; 1:tapering is enabled on point sites
    end
    
    methods
        
        function obj = stem_model(stem_data,stem_par)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %stem_data - [stem_data object]  (1x1)
            %stem_par  - [stem_par object]   (1x1)
            %
            %OUTPUT
            %obj       - [stem_model object] (1x1)
            if nargin<2
                error('Not enough input arguments');
            end
            obj.stem_data=stem_data;
            obj.stem_par=stem_par;
            if not(length(obj.stem_par.beta)==size(obj.stem_data.X_beta,2))
                error(['The length of beta in stem_par must be equal to ',num2str(size(obj.stem_data.X_beta,2))]);
            end
            
            if not(isempty(stem_data.stem_gridlist_b))
                if not(isempty(stem_data.stem_gridlist_b.tap))
                    obj.tapering_b=1;
                else
                    obj.tapering_b=0;
                end
            end
            if not(isempty(stem_data.stem_gridlist_p.tap))
                obj.tapering_p=1;
            else
                obj.tapering_p=0;
            end
            if not(isempty(stem_data.stem_gridlist_b))
                obj.tapering=obj.tapering_p|obj.tapering_b;
            else
                obj.tapering=obj.tapering_p;
            end
        end
        
        function [aj_bp,aj_p,aj_z] = get_aj(obj)
            %DESCRIPTION: provides the vector aj_bp, aj_p and aj_z used in the EM estimation
            %
            %INPUT
            %obj   - [stem_model object] (1x1)
            %
            %OUTPUT
            %aj_bp - [double]            (Nx1) is the diagonal of the NxN diagonal matrix alpha_bp*J_bp
            %aj_p  - [double]            (NxK) are the diagonals of the NxN diagonal matrices alpha_p*J_p 
            %aj_z  - [double]            (Nx1) is the diagonal of the NxN diagonal matrix alpha_z*J_z
            %
            %NOTE
            %The elements of aj_p and aj_z from Np+1 to N are all zeros. This allows the
            %use of stem_misc.D_apply both for the pixel data and the point
            %level data avoiding the use of J_bp, J_p and J_z
            if not(isempty(obj.stem_data.stem_varset_b))
                aj_bp=zeros(obj.stem_data.N,1);
                blocks=[0 cumsum(obj.stem_data.dim)];
                for i=1:obj.stem_data.nvar
                    aj_bp(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
            else
                aj_bp=[];
            end
            if obj.stem_par.k>0
                aj_p=zeros(obj.stem_data.N,obj.stem_par.k);
                blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                for k=1:obj.stem_par.k
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        aj_p(blocks(i)+1:blocks(i+1),k)=obj.stem_par.alpha_p(i,k);    
                    end
                end
            else
                aj_p=[];
            end
            if obj.stem_data.model_type==1
                dim=obj.stem_data.dim;
                if obj.stem_data.model_subtype==1
                    aj_z=kron(obj.stem_par.alpha_z,ones(dim(1),1));
                else
                    aj_z=[];         
                    for i=1:obj.stem_par.p
                        aj_z=cat(1,aj_z,obj.stem_par.alpha_z(i)*ones(dim(i),1));
                    end
                end
                z=zeros(obj.stem_data.N-length(aj_z),1);
                aj_z=[aj_z;z];
            else
                aj_z=[];
            end
        end
        
        function [aj_bp_b,j_b] = get_jbp(obj,r)
            %DESCRIPTION: provides the vectors aj_bp_b and j_b used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %
            %OUTPUT
            %aj_bp_b - [double]            (Nx1) is the vector with elements equal to alpha_bp(r) only for the sites of the r-th variable
            %j_b     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            aj_bp_b=zeros(obj.stem_data.N,1);
            j_b=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.nvar
                if i==r
                    j_b(blocks(i)+1:blocks(i+1))=1;
                    aj_bp_b(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_bp(i);
                end
            end
        end
        
        function [aj_p_bs,j_p] = get_jp(obj,r,s)
            %DESCRIPTION: provides the vectors aj_p_bs and j_p used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and q
            %s       - [integer]           (1x1) is the index between 1 and K
            %
            %OUTPUT
            %aj_p_bs - [double]            (Nx1) is the vector with elements equal to alpha_p(i,s) only for the sites of the r-th variable
            %j_p     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            aj_p_bs=zeros(obj.stem_data.N,1);
            j_p=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_data.stem_varset_p.nvar
                if i==r
                    j_p(blocks(i)+1:blocks(i+1))=1;
                    aj_p_bs(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_p(i,s); 
                end
            end
        end  
        
        function [aj_z_p,j_z] = get_jz(obj,r)
            %DESCRIPTION: provides the vectors aj_z_p and j_z used in the EM estimation
            %
            %INPUT
            %obj     - [stem_model object] (1x1)
            %r       - [integer]           (1x1) is the index between 1 and p
            %
            %OUTPUT
            %aj_z_p  - [double]            (Nx1) is the vector with elements equal to alpha_z(r) only for the sites of the r-th variable
            %j_z     - [double]            (Nx1) is the vector with elements equal to 1 only for the sites of the r-th variable
            
            aj_z_p=zeros(obj.stem_data.N,1);
            j_z=zeros(obj.stem_data.N,1);
            blocks=[0 cumsum(obj.stem_data.dim)];
            for i=1:obj.stem_par.p
                if i==r
                    j_z(blocks(i)+1:blocks(i+1))=1;
                    aj_z_p(blocks(i)+1:blocks(i+1))=obj.stem_par.alpha_z(i);
                end
            end
        end        
        
        function [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,sigma_Z,sigma_eta,G_tilde_diag,aj_bp,aj_p,aj_z,M] = get_sigma(obj,sigma_W_b)
            %DESCRIPTION: provides the variance-covariance matrices and some vectors that are used in the EM algorithm
            %
            %INPUT
            %
            %obj         - [stem_model object] (1x1)
            %<sigma_W_b> - [double]            (NbxNb) (default: []) variance-covariance matrix of W_b. It is provided as input argument during kriging when sigma_W_b does not change across blocks
            %
            %OUTPUT
            %
            %sigma_eps   - [double]            (NxN) variance-covariance matrix of epsilon
            %sigma_W_b   - [double]            (N_bxN_b) variance-covariance matrix of W_b
            %sigma_W_p   - [double]            {K}(N_pxN_p) variance-covariance matrices of the K W_p_i
            %sigma_geo   - [double]            (NxN) variance-covariance matrix of the sum of all the geostatistical components (Z excluded and epsilon included)
            %sigma_Z     - [double]            (pxp) variance-covariance of Z
            %sigma_eta   - [double]            (r x r) variance-covariance matrix of eta when model_type=1
            %G_tilde_diag- [double]            (r x 1) diagonal of the G_tilde matrix when model_type=1
            %aj_bp       - [double]            (Nx1) see the details of the method get_jbp
            %aj_p        - [double]            (Nx1) see the details of the method get_jp
            %aj_z        - [double]            (Nx1) see the details of the method get_jz
            %M           - [integer >0]        (N_px1) see the details of the method update_M of the class stem_data
            %
            %NOTE
            %sigma_geo is provided only if it is time-invariant otherwise it is evaluated at each step of the EM algorithm
            disp('    Marginal variance-covariance matrices evaluation started...');
            ct1=clock;
            if nargin==1
                sigma_W_b=[];
            end
            
            nvar=obj.stem_data.nvar;
            N=obj.stem_data.N;
            M=obj.stem_data.M;
            dim=obj.stem_data.dim;
            
            %sigma_eps
            d=[];
            for i=1:nvar
                d=cat(1,d,repmat(obj.stem_par.sigma_eps(i,i),dim(i),1));
            end
            
            I=1:length(d);
            sigma_eps=sparse(I,I,d);
            
            %sigma_W_b
            if not(isempty(obj.stem_data.stem_varset_b))
                if not(isempty(obj.stem_data.stem_varset_b.X_bp))
                    if (nargin==1)||isempty(sigma_W_b)
                        if not(obj.tapering_b)
                            sigma_W_b=zeros(obj.stem_data.stem_varset_b.N);
                        end
                        blocks=[0 cumsum(obj.stem_data.stem_varset_b.dim)];
                        if obj.stem_par.pixel_correlated
                            if obj.tapering_b
                                I=zeros(nnz(obj.stem_data.DistMat_b),1);
                                J=zeros(nnz(obj.stem_data.DistMat_b),1);
                                elements=zeros(nnz(obj.stem_data.DistMat_b),1);
                            end
                            idx=0;
                            for j=1:obj.stem_data.stem_varset_b.nvar
                                for i=j:obj.stem_data.stem_varset_b.nvar
                                    idx_r=blocks(i)+1:blocks(i+1);
                                    idx_c=blocks(j)+1:blocks(j+1);
                                    if not(obj.tapering_b)
                                        sigma_W_b(idx_r,idx_c)=obj.stem_par.v_b(i,j)*stem_misc.correlation_function(...
                                            obj.stem_par.theta_b,obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_par.correlation_type);
                                        if not(i==j)
                                            sigma_W_b(idx_c,idx_r)=sigma_W_b(idx_r,idx_c)';
                                        end
                                    else
                                        corr_result=stem_misc.correlation_function(obj.stem_par.theta_b,obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_par.correlation_type);
                                        weights=stem_misc.wendland(obj.stem_data.DistMat_b(idx_r,idx_c),obj.stem_data.stem_gridlist_b.tap);
                                        corr_result.correlation=obj.stem_par.v_b(i,j)*corr_result.correlation.*weights;
                                        siz=length(corr_result.I);
                                        I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                        J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                        elements(idx+1:idx+siz)=corr_result.correlation;
                                        idx=idx+siz;
                                        if not(i==j)
                                            I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                            J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                            elements(idx+1:idx+siz)=corr_result.correlation;
                                            idx=idx+siz;
                                        end
                                    end
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        else
                            if obj.tapering_b
                                idx=0;
                                nonzeros=0;
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    nonzeros=nonzeros+nnz(obj.stem_data.DistMat_b(blocks(i)+1:blocks(i+1),blocks(i)+1:blocks(i+1)));
                                end
                                I=zeros(nonzeros,1);
                                elements=zeros(nonzeros,1);
                            end
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                idx_rc=blocks(i)+1:blocks(i+1);
                                if not(obj.tapering_b)
                                    sigma_W_b(idx_rc,idx_rc)=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                else
                                    corr_result=stem_misc.correlation_function(obj.stem_par.theta_b(:,i),obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_par.correlation_type);
                                    weights=stem_misc.wendland(obj.stem_data.DistMat_b(idx_rc,idx_rc),obj.stem_data.stem_gridlist_b.tap);
                                    corr_result.correlation=obj.stem_par.v_b(i,i)*corr_result.correlation.*weights;
                                    siz=length(corr_result.I);
                                    I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    J(idx+1:idx+siz)=corr_result.J+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                            if obj.tapering_b
                                sigma_W_b=sparse(I,J,elements);
                            end
                        end
                    end
                end
            end
            clear I
            clear J
            clear elements
            clear weights
            clear corr_result
            
            %sigma_W_p
            if obj.stem_par.k>0
                blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                sigma_W_p=cell(obj.stem_par.k,1);
                for k=1:obj.stem_par.k
                    if obj.tapering_p
                        I=zeros(nnz(obj.stem_data.DistMat_p),1);
                        J=zeros(nnz(obj.stem_data.DistMat_p),1);
                        elements=zeros(nnz(obj.stem_data.DistMat_p),1);
                    else
                        sigma_W_p{k}=obj.stem_data.DistMat_p;
                    end
                    idx=0;
                    for j=1:obj.stem_data.stem_varset_p.nvar
                        for i=j:obj.stem_data.stem_varset_p.nvar
                            idx_r=blocks(i)+1:blocks(i+1);
                            idx_c=blocks(j)+1:blocks(j+1);
                            if not(obj.tapering_p)
                                sigma_W_p{k}(idx_r,idx_c)=obj.stem_par.v_p(i,j,k)*stem_misc.correlation_function(...
                                    obj.stem_par.theta_p(:,k),obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                                if not(i==j)
                                    sigma_W_p{k}(idx_c,idx_r)=sigma_W_p{k}(idx_r,idx_c)';
                                end
                            else
                                corr_result=stem_misc.correlation_function(obj.stem_par.theta_p(:,k),obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                                weights=stem_misc.wendland(obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_data.stem_gridlist_p.tap);
                                corr_result.correlation=obj.stem_par.v_p(i,j,k)*corr_result.correlation.*weights;
                                siz=length(corr_result.I);
                                I(idx+1:idx+siz)=corr_result.I+blocks(i);
                                J(idx+1:idx+siz)=corr_result.J+blocks(j);
                                elements(idx+1:idx+siz)=corr_result.correlation;
                                idx=idx+siz;
                                if not(i==j)
                                    I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                    J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                    elements(idx+1:idx+siz)=corr_result.correlation;
                                    idx=idx+siz;
                                end
                            end
                        end
                    end
                    if obj.tapering_p
                        sigma_W_p{k}=sparse(I,J,elements);
                    end
                end
            else
                sigma_W_p=[];
            end
            clear I
            clear J
            clear elements
            clear weights
            clear corr_result
            
            %sigma_eta
            if obj.stem_data.model_type==1
                if not(obj.tapering_p)
                    sigma_eta=zeros(size(obj.stem_par.G));
                end
                
                if obj.stem_data.model_subtype==0
                    blocks=[0 cumsum(obj.stem_data.stem_varset_p.dim)];
                else
                    blocks=[0 cumsum(repmat(obj.stem_data.stem_varset_p.dim(1),[1,obj.stem_par.p]))];
                end
                
                if obj.tapering_p
                    if obj.stem_data.model_subtype==0
                        I=zeros(nnz(obj.stem_data.DistMat_p),1);
                        J=zeros(nnz(obj.stem_data.DistMat_p),1);
                        elements=zeros(nnz(obj.stem_data.DistMat_p),1);
                    else
                        I=zeros(nnz(obj.stem_data.DistMat_z),1);
                        J=zeros(nnz(obj.stem_data.DistMat_z),1);
                        elements=zeros(nnz(obj.stem_data.DistMat_z),1);
                    end
                end
                
                idx=0;
                for j=1:obj.stem_par.p
                    for i=j:obj.stem_par.p
                        idx_r=blocks(i)+1:blocks(i+1);
                        idx_c=blocks(j)+1:blocks(j+1);
                        if not(obj.tapering_p)
                            if obj.stem_data.model_subtype==0
                                sigma_eta(idx_r,idx_c)=obj.stem_par.v_z(i,j)*stem_misc.correlation_function(...
                                    obj.stem_par.theta_z,obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                            else
                                sigma_eta(idx_r,idx_c)=obj.stem_par.v_z(i,j)*stem_misc.correlation_function(...
                                    obj.stem_par.theta_z,obj.stem_data.DistMat_z(idx_r,idx_c),obj.stem_par.correlation_type);
                            end
                            if not(i==j)
                                sigma_eta(idx_c,idx_r)=sigma_eta(idx_r,idx_c)';
                            end
                        else
                            if obj.stem_data.model_subtype==0
                                corr_result=stem_misc.correlation_function(obj.stem_par.theta_z,obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_par.correlation_type);
                                weights=stem_misc.wendland(obj.stem_data.DistMat_p(idx_r,idx_c),obj.stem_data.stem_gridlist_p.tap);
                            else
                                corr_result=stem_misc.correlation_function(obj.stem_par.theta_z,obj.stem_data.DistMat_z(idx_r,idx_c),obj.stem_par.correlation_type);
                                weights=stem_misc.wendland(obj.stem_data.DistMat_z(idx_r,idx_c),obj.stem_data.stem_gridlist_p.tap);
                            end
                            corr_result.correlation=obj.stem_par.v_z(i,j)*corr_result.correlation.*weights;
                            siz=length(corr_result.I);
                            I(idx+1:idx+siz)=corr_result.I+blocks(i);
                            J(idx+1:idx+siz)=corr_result.J+blocks(j);
                            elements(idx+1:idx+siz)=corr_result.correlation;
                            idx=idx+siz;
                            if not(i==j)
                                I(idx+1:idx+siz)=corr_result.J+blocks(j);
                                J(idx+1:idx+siz)=corr_result.I+blocks(i);
                                elements(idx+1:idx+siz)=corr_result.correlation;
                                idx=idx+siz;
                            end
                        end
                    end
                end
                if obj.tapering_p
                    sigma_eta=sparse(I,J,elements);
                end
            else
                sigma_eta=[];
            end
            
            %sigma_geo
            sigma_geo=[];
            [aj_bp,aj_p,aj_z]=obj.get_aj;
            if not(obj.stem_data.X_tv)
                %time invariant case
                if not(isempty(obj.stem_data.stem_varset_b))
                    sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),obj.stem_data.X_bp(:,1,1),'b'),aj_bp,'b');
                end
                if obj.stem_par.k>0
                    if isempty(sigma_geo)
                        %se manca il remoto allora sigma_geo non � stata
                        %ancora allocata
                        if obj.tapering
                            sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                        else
                            sigma_geo=zeros(N);
                        end
                    end
                    for k=1:obj.stem_par.k
                        sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{k},obj.stem_data.X_p(:,1,1,k),'b'),aj_p(:,k),'b');
                    end
                end
                if isempty(sigma_geo)
                    sigma_geo=sigma_eps;
                else
                    sigma_geo=sigma_geo+sigma_eps;
                end
            end
            
            if obj.stem_par.p>0
                if not(obj.stem_data.model_type==1)
                    sigma_Z=(eye(obj.stem_par.p^2)-kron(obj.stem_par.G,obj.stem_par.G))\obj.stem_par.sigma_eta(:);%variance of a VAR (Lutkepohl pag.27)
                    sigma_Z=reshape(sigma_Z,obj.stem_par.p,obj.stem_par.p);
                    G_tilde_diag=[];
                else
                    dim=obj.stem_data.dim;
                    if obj.stem_data.model_subtype==1
                        G_tilde_diag=kron(diag(obj.stem_par.G),ones(dim(1),1));
                    else
                        G_tilde_diag=[];
                        for i=1:obj.stem_par.p
                            G_tilde_diag=cat(1,G_tilde_diag,obj.stem_par.G(i,i)*ones(dim(i),1));
                        end
                    end
                    sigma_Z=reshape(1./(1-kron(G_tilde_diag,G_tilde_diag)).*sigma_eta(:),length(G_tilde_diag),length(G_tilde_diag));
                end
            else
                sigma_Z=[];
                G_tilde_diag=[];
            end
            ct2=clock;
            disp(['    Marginal variance-covariance matrices evaluation ended in ',stem_misc.decode_time(etime(ct2,ct1))]);
        end
        
        function simulate(obj,nan_rate,nan_pattern_par)
            %DESCRIPTION: is the front-end method to simulate data from this stem_model object
            %
            %INPUT
            %
            %obj               - [stem_model object] (1x1)
            %<nan_rate>        - [double [0,1]]      (Sx1) (default: []) missing rate of the simulated data
            %<nan_pattern_par> - [double >0]         (Sx1) (default: []) (UoM: km) parameter of the exponential spatial correlation function used to define the spatial pattern of the missing data
            %
            %OUTPUT
            %
            %none: an object of class stem_sim is created. The method stem_sim.simulate actually simulates the data and it fills the stem_data object
            st_sim=stem_sim(obj);
            if nargin<2
                nan_rate=[];
                nan_pattern_par=[];
            end
            if nargin==2
                error('nan_pattern_par must also be provided');
            end
            st_sim.simulate(nan_rate,nan_pattern_par);
            obj.stem_data.simulated=1;
        end        
        
        function EM_estimate(obj,stem_EM_options)
            %DESCRIPTION: front-end method for EM estimation of this model
            %
            %INPUT
            %
            %obj                - [stem_model object]      (1x1)
            %stem_EM_options    - [stem_EM_options object] (1x)
            %
            %OUTPUT
            %
            %none: the stem_EM object is created. The method stem_EM.estimate is used to estimate the model and it updates the stem_par object

            obj.set_system_size;
            
            standardized=1;
            if not(obj.stem_data.stem_varset_p.standardized)
                standardized=0;
            end
            if not(isempty(obj.stem_data.stem_varset_b))
                if not(obj.stem_data.stem_varset_b.standardized)
                    standardized=0;
                end
            end
            if not(standardized) && not(obj.stem_data.model_type==2||obj.stem_data.model_type==3)
                warning('All X and Y are standardized in order to avoid numerical stability problems. Consider to explicitly use the method standardize of class stem_data. The estimated model parameters will reflect standardization.');
                obj.stem_data.standardize;
                obj.stem_par.beta = obj.get_beta0();
                obj.stem_par_initial.beta=obj.stem_par.beta;
            end
            if not(isempty(obj.stem_data.stem_crossval))
                disp('Data modification for cross-validation started...');
                for i=1:length(obj.stem_data.stem_crossval.variable_name)
                    obj.stem_data.stem_crossval.stem_crossval_result{i}=stem_crossval_result();
                    
                    idx_var=obj.stem_data.stem_varset_p.get_Y_index(obj.stem_data.stem_crossval.variable_name{i});
                    if isempty(idx_var)
                        error(['Cross-validation variable',obj.stem_data.stem_crossval.variable_name{i},'not found.']);
                    end
                    
                    %recover the indices of the cross-validation sites
                    indices=obj.stem_data.stem_crossval.indices{i};
                    Y={obj.stem_data.stem_varset_p.Y{idx_var}(indices,:)};
                    Y_name=obj.stem_data.stem_varset_p.Y_name(idx_var);
                    if not(isempty(obj.stem_data.stem_varset_p.X_bp))
                        X_bp={obj.stem_data.stem_varset_p.X_bp{idx_var}(indices,:,:)};
                        X_bp_name=obj.stem_data.stem_varset_p.X_bp(idx_var);
                    else
                        X_bp={};
                        X_bp_name={};
                    end
                    if not(isempty(obj.stem_data.stem_varset_p.X_beta))
                        X_beta={obj.stem_data.stem_varset_p.X_beta{idx_var}(indices,:,:)};
                        X_beta_name=obj.stem_data.stem_varset_p.X_beta_name(idx_var);
                    else
                        X_beta={};
                        X_beta_name={};
                    end
                    if not(isempty(obj.stem_data.stem_varset_p.X_z))
                        X_z={obj.stem_data.stem_varset_p.X_z{idx_var}(indices,:,:)};
                        X_z_name=obj.stem_data.stem_varset_p.X_z_name(idx_var);
                    else
                        X_z={};
                        X_z_name={};
                    end
                    if not(isempty(obj.stem_data.stem_varset_p.X_p))
                        X_p={obj.stem_data.stem_varset_p.X_p{idx_var}(indices,:,:,:)};
                        X_p_name=obj.stem_data.stem_varset_p.X_p_name(idx_var);
                    else
                        X_p={};
                        X_p_name={};
                    end
                    
                    %set the cross_mindistance vector
                    if not(isempty(obj.stem_data.DistMat_p))
                        dim=obj.stem_data.dim;
                        blocks=[0 cumsum(dim)];
                        temp_dist=obj.stem_data.DistMat_p(blocks(idx_var)+1:blocks(idx_var+1),blocks(idx_var)+1:blocks(idx_var+1));
                        temp_dist=temp_dist(indices,:);
                        temp_dist(:,indices)=[];
                        obj.stem_data.stem_crossval.stem_crossval_result{i}.min_distance=min(temp_dist,[],2);
                        clear temp_dist
                    end
                    
                    obj.stem_data.stem_crossval.stem_varset{i}=stem_varset(Y,Y_name,X_bp,X_bp_name,X_beta,X_beta_name,X_z,X_z_name,X_p,X_p_name);
                    obj.stem_data.stem_crossval.stem_gridlist{i}=stem_gridlist();
                    coordinate=obj.stem_data.stem_gridlist_p.grid{idx_var}.coordinate;
                    st_grid=stem_grid(coordinate(indices,:),'deg','sparse','point');
                    obj.stem_data.stem_crossval.stem_gridlist{i}.add(st_grid);
                    %remove the cross-validation data from the estimation dataset
                    obj.stem_data.site_crop(obj.stem_data.stem_crossval.type{i},obj.stem_data.stem_crossval.variable_name{i},indices,1);
                end
                obj.cross_validation=1;
                disp('Data modification ended.');
            else
                obj.cross_validation=0;
            end
            
            st_EM=stem_EM(obj,stem_EM_options);
            %set the current parameter value with the estimated initial value
            obj.stem_par=obj.stem_par_initial;
            if isempty(stem_EM_options.path_distributed_computing)
                obj.stem_EM_result=st_EM.estimate();
            else
                obj.stem_EM_result=st_EM.estimate_parallel(stem_EM_options.path_distributed_computing);
            end
            obj.estimated=1;
            if obj.cross_validation
                obj.fill_crosval_result();
            end
        end
          
        function fill_crosval_result(obj)
            if obj.cross_validation
                for i=1:length(obj.stem_data.stem_crossval.variable_name)
                    disp(['Kriging over cross-validation sites of variable ',obj.stem_data.stem_crossval.variable_name{i}]);
                    
                    idx_var=obj.stem_data.stem_varset_p.get_Y_index(obj.stem_data.stem_crossval.variable_name{i});
                    if isempty(idx_var)
                        error(['Cross-validation variable',obj.stem_data.stem_crossval.variable_name{i},'not found.']);
                    end
                    
                    st_krig=stem_krig(obj);
                    block_size=1000;
                    back_transform=0;
                    no_varcov=0;
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.stem_krig_result=st_krig.kriging(obj.stem_data.stem_crossval.variable_name{i},obj.stem_data.stem_crossval.stem_gridlist{i}.grid{1},block_size,[],[],back_transform,no_varcov,i);
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.res=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}-obj.stem_data.stem_crossval.stem_crossval_result{i}.stem_krig_result.y_hat;
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.mse=nanvar(obj.stem_data.stem_crossval.stem_crossval_result{i}.res');
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.mse_time=nanvar(obj.stem_data.stem_crossval.stem_crossval_result{i}.res);
                    
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.relative_mse=obj.stem_data.stem_crossval.stem_crossval_result{i}.mse./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}');
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.relative_mse_time=obj.stem_data.stem_crossval.stem_crossval_result{i}.mse_time./nanvar(obj.stem_data.stem_crossval.stem_varset{i}.Y{1});
                    
                    if (obj.stem_data.stem_varset_p.standardized)&&not(obj.stem_data.stem_varset_p.log_transformed)
                        s=obj.stem_data.stem_varset_p.Y_stds{idx_var};
                        m=obj.stem_data.stem_varset_p.Y_means{idx_var};
                        y_hat_back=obj.stem_data.stem_crossval.stem_crossval_result{i}.stem_krig_result.y_hat*s+m;
                        y=obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m;
                    end
                    
                    if (obj.stem_data.stem_varset_p.standardized)&&(obj.stem_data.stem_varset_p.log_transformed)
                        s=obj.stem_data.stem_varset_p.Y_stds{idx_var};
                        m=obj.stem_data.stem_varset_p.Y_means{idx_var};
                        y_hat_back=obj.stem_data.stem_crossval.stem_crossval_result{i}.stem_krig_result.y_hat;
                        var_y_hat=obj.stem_data.stem_crossval.stem_crossval_result{i}.stem_krig_result.diag_Var_y_hat;
                        y_hat_back=exp(y_hat_back*s+m+(var_y_hat*s^2)/2);
                        y=exp(obj.stem_data.stem_crossval.stem_varset{i}.Y{1}*s+m);
                    end
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.res_back=y-y_hat_back;
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.y_back=y;
                    obj.stem_data.stem_crossval.stem_crossval_result{i}.y_hat_back=y_hat_back;
                end
            else
                disp('The stem_model object does not include cross validation information');
            end
        end
          
        function print(obj)
            %DESCRIPTION: print the information on the estimation result
            %
            %INPUT
            %obj  - [stem_model object]   (1x1) the stem_data object
            %
            %OUTPUT
            %
            %none: the information is printed in the command window          
            
            stem_misc.disp_star('Model estimation results')   
            if obj.estimated
                if obj.tapering
                    disp('* Tapering is enabled.');    
                    if not(isempty(obj.tapering_p))
                        disp(['  Point data tapering: ',num2str(obj.stem_data.stem_gridlist_p.tap),' km']);
                    else
                        disp('   Tapering is NOT enabled on point data');
                    end
                    if not(isempty(obj.tapering_b))
                        disp(['  Pixel data tapering: ',num2str(obj.stem_data.stem_gridlist_b.tap),' km']);
                    else
                        disp('  Tapering is NOT enabled on pixel data');
                    end
                else
                    disp('* Tapering is not enabled');
                end
                disp(' ');
                if not(isempty(obj.stem_EM_result.logL))
                    disp(['* Observed data log-likelihood: ',num2str(obj.stem_EM_result.logL,'%05.3f')]);
                else
                    disp('* Observed data log-likelihood: not computed. Use the method set_logL of the class stem_model.');
                end
                disp(' ');
                counter=1;
                if not(isempty(obj.stem_par.beta))
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        disp(['* Beta coefficients related to the point variable ',obj.stem_data.stem_varset_p.Y_name{i}]);
                        output=cell(length(obj.stem_data.stem_varset_p.X_beta_name{i})+1,3);
                        output{1,1}='Loading coefficient';
                        output{1,2}='Value';
                        output{1,3}='Std';
                        for j=1:length(obj.stem_data.stem_varset_p.X_beta_name{i})
                            output{j+1,1}=obj.stem_data.stem_varset_p.X_beta_name{i}{j};
                            output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{j+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                    end
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_beta))
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                disp(['* Beta coefficients related to the pixel variable ',obj.stem_data.stem_varset_b.Y_name{i}]);
                                output=cell(length(obj.stem_data.stem_varset_b.X_beta_name{i})+1,3);
                                output{1,1}='Loading coefficient';
                                output{1,2}='Value';
                                output{1,3}='Std';
                                for j=1:length(obj.stem_data.stem_varset_b.X_beta_name{i})
                                    output{j+1,1}=obj.stem_data.stem_varset_b.X_beta_name{i}{j};
                                    output{j+1,2}=num2str(obj.stem_par.beta(counter),'%+05.3f');
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{j+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                                    else
                                        output{j+1,3}='Not computed';
                                    end
                                    counter=counter+1;
                                end
                                disp(output);
                            end
                        end
                    end
                end
                output=cell(obj.stem_data.stem_varset_p.nvar+1,3);
                disp('* Sigma_eps diagonal elements (Variance)')
                output{1,1}='Variable';
                output{1,2}='Value';
                output{1,3}='Std';
                for i=1:obj.stem_data.stem_varset_p.nvar
                    output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                    output{i+1,2}=num2str(obj.stem_par.sigma_eps(i,i),'%05.3f');
                    if not(isempty(obj.stem_EM_result.varcov))
                        output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                    else
                        output{i+1,3}='Not computed';
                    end
                    counter=counter+1;
                end
                if not(isempty(obj.stem_data.stem_varset_b))
                    delta=obj.stem_data.stem_varset_p.nvar;
                    for i=1:obj.stem_data.stem_varset_b.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.sigma_eps(i+delta,i+delta),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                end
                disp(output);
                if not(isempty(obj.stem_data.stem_varset_b))
                    disp('* alpha_bp elements')
                    output=cell(obj.stem_data.stem_varset_p.nvar+1,3);
                    output{1,1}='Variable';
                    output{1,2}='Value';
                    output{1,3}='Std';
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        output{i+1,2}=num2str(obj.stem_par.alpha_bp(i),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    delta=obj.stem_data.stem_varset_p.nvar;
                    for i=1:obj.stem_data.stem_varset_b.nvar
                        output{i+1+delta,1}=obj.stem_data.stem_varset_b.Y_name{i};
                        output{i+1+delta,2}=num2str(obj.stem_par.alpha_bp(i+delta),'%+05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{i+1+delta,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{i+1+delta,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    output=[];
                    if obj.stem_par.pixel_correlated
                        disp('* Pixel data are cross-correlated.');
                        disp(' ');
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        output{2,1}='Theta_b';
                        output{2,2}=num2str(obj.stem_par.theta_b(1),'%05.3f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{2,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                        else
                            output{2,3}='Not computed';
                        end
                        counter=counter+1;
                        disp(output);
                        output=cell(obj.stem_data.stem_varset_b.nvar+1,obj.stem_data.stem_varset_b.nvar+1);
                        disp('* v_b matrix:');
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            for j=i+1:obj.stem_data.stem_varset_b.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_b(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_b(i,j),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end
                        disp(output);
                    else
                        disp('* Pixel data are NOT cross-correlated.');
                        disp(' ');
                        disp('* Theta_b elements:');
                        output=cell(obj.stem_data.stem_varset_b.nvar+1,3);
                        output{1,1}='Variable';
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        for i=1:obj.stem_data.stem_varset_b.nvar
                            output{i+1,1}=obj.stem_data.stem_varset_b.Y_name{i};
                            output{i+1,2}=num2str(obj.stem_par.theta_b(i),'%05.3f');
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f');
                            else
                                output{i+1,3}='Not computed';
                            end
                            counter=counter+1;
                        end
                        disp(output);
                    end
                end
                output=cell(obj.stem_data.stem_varset_p.nvar*2,obj.stem_par.k+1);
                if obj.stem_par.k>0
                    disp(['* ',num2str(obj.stem_par.k),' fine-scale coregionalization components w_p']);
                    disp(' ');
                    disp('* alpha_p elements:')
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        output{i*2,1}=obj.stem_data.stem_varset_p.Y_name{i};
                        for k=1:obj.stem_par.k
                            output{i*2-1,k+1}=obj.stem_data.stem_varset_p.X_p_name{i}{k};
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i*2,k+1}=[num2str(obj.stem_par.alpha_p(i,k),'%+05.3f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f'),')'];
                            else
                                output{i*2,k+1}=num2str(obj.stem_par.alpha_p(i,k),'%+05.3f');
                            end
                            counter=counter+1;
                        end
                    end
                    disp(output);
                    disp('* theta_p elements:');
                    output=cell(1,3);
                    output{1,1}='Coreg. component';
                    output{1,2}='Value [km]';
                    output{1,3}='Std [km]';
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end
                        output{k+1,1}=[num2str(k),postfix];
                        output{k+1,2}=num2str(obj.stem_par.theta_p(:,k),'%06.2f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{k+1,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f');
                        else
                            output{k+1,3}='Not computed';
                        end
                        counter=counter+1;
                    end
                    disp(output);
                    for k=1:obj.stem_par.k
                        if k==1
                            postfix='st';
                        end
                        if k==2
                            postfix='nd';
                        end
                        if k==3
                            postfix='rd';
                        end
                        if k>3
                            postfix='th';
                        end                        
                        disp(['* v_p matrix for the ',num2str(k),postfix,' coreg. component:']);
                        output=cell(obj.stem_data.stem_varset_p.nvar+1,obj.stem_data.stem_varset_p.nvar+1);
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            output{1,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=i+1:obj.stem_data.stem_varset_p.nvar
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_p(i,j,k),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_p(i,j,k),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end                        
                        disp(output);
                    end
                end
                if obj.stem_par.p>0
                    output=cell(obj.stem_par.p+1,obj.stem_par.p+1);
                    disp('* Transition matrix G:');
                    c=1;
                    for i=1:obj.stem_data.stem_varset_p.nvar
                        for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                            output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                            output{c+1,1}=output{1,c+1};
                            c=c+1;
                        end
                    end
                    if not(isempty(obj.stem_data.stem_varset_b))
                        if not(isempty(obj.stem_data.stem_varset_b.X_z))
                            for i=1:obj.stem_data.stem_varset_b.nvar
                                for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                    output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                    output{c+1,1}=output{1,c+1};
                                    c=c+1;
                                end
                            end
                        end
                    end
                    if (obj.stem_par.time_diagonal) || (obj.stem_data.model_type==1)
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if i==j
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,i+1}=[num2str(obj.stem_par.G(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,i+1}=num2str(obj.stem_par.G(i,i),'%+05.2f');
                                    end
                                    counter=counter+1;
                                else
                                    output{i+1,j+1}='0';
                                end
                            end
                        end
                    else
                        for j=1:obj.stem_par.p
                            for i=1:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.G(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.G(i,j),'%+05.2f');
                                end
                                counter=counter+1;
                            end
                        end
                    end
                    disp(output);
                    if not(obj.stem_data.model_type==1)
                        disp('* Sigma_eta matrix:');
                        c=1;
                        output=cell(obj.stem_par.p+1,obj.stem_par.p);
                        for i=1:obj.stem_data.stem_varset_p.nvar
                            for j=1:size(obj.stem_data.stem_varset_p.X_z{i},2)
                                output{1,c+1}=[obj.stem_data.stem_varset_p.Y_name{i},' - ',obj.stem_data.stem_varset_p.X_z_name{i}{j}];
                                output{c+1,1}=output{1,c+1};
                                c=c+1;
                            end
                        end
                        if not(isempty(obj.stem_data.stem_varset_b))
                            if not(isempty(obj.stem_data.stem_varset_b.X_z))
                                for i=1:obj.stem_data.stem_varset_b.nvar
                                    for j=1:size(obj.stem_data.stem_varset_b.X_z{i},2)
                                        output{1,c+1}=[obj.stem_data.stem_varset_b.Y_name{i},' - ',obj.stem_data.stem_varset_b.X_z_name{i}{j}];
                                        output{c+1,1}=output{1,c+1};
                                        c=c+1;
                                    end
                                end
                            end
                        end
                        if obj.stem_par.time_diagonal
                            for j=1:obj.stem_par.p
                                for i=1:obj.stem_par.p
                                    if i==j
                                        if not(isempty(obj.stem_EM_result.varcov))
                                            output{i+1,i+1}=[num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                        else
                                            output{i+1,i+1}=num2str(obj.stem_par.sigma_eta(i,i),'%+05.2f');
                                        end
                                        counter=counter+1;
                                    else
                                        output{i+1,j+1}='0';
                                    end
                                end
                            end
                        else
                            for j=1:obj.stem_par.p
                                for i=j:obj.stem_par.p
                                    if not(isempty(obj.stem_EM_result.varcov))
                                        output{i+1,j+1}=[num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    else
                                        output{i+1,j+1}=num2str(obj.stem_par.sigma_eta(i,j),'%+05.2f');
                                    end
                                    output{j+1,i+1}=output{i+1,j+1};
                                    counter=counter+1;
                                end
                            end
                        end
                        disp(output);
                    else
                        disp('* Fine-scale coregionalization components z');
                        disp(' ');
                        disp('* alpha_z elements:')
                        output=cell(obj.stem_par.p*2,2);
                        for i=1:obj.stem_par.p
                            output{i*2,1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i*2-1,2}=cell2mat(obj.stem_data.stem_varset_p.X_z_name{i});
                            if not(isempty(obj.stem_EM_result.varcov))
                                output{i*2,2}=[num2str(obj.stem_par.alpha_z(i),'%+05.3f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.3f'),')'];
                            else
                                output{i*2,2}=num2str(obj.stem_par.alpha_z(i),'%+05.3f');
                            end
                            counter=counter+1;
                        end
                        disp(output);
                        disp('* theta_z elements:');
                        output=cell(1,3);
                        output{1,1}='Coreg. component';
                        output{1,2}='Value [km]';
                        output{1,3}='Std [km]';
                        
                        output{2,2}=num2str(obj.stem_par.theta_z,'%06.2f');
                        if not(isempty(obj.stem_EM_result.varcov))
                            output{2,3}=num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%05.2f');
                        else
                            output{2,3}='Not computed';
                        end
                        counter=counter+1;
                        
                        disp(output);
                        
                        disp('* v_z matrix for the coreg. component:');
                        output=cell(obj.stem_par.p+1,obj.stem_par.p+1);
                        for i=1:obj.stem_par.p
                            output{1,i+1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,1}=obj.stem_data.stem_varset_p.Y_name{i};
                            output{i+1,i+1}=num2str(1,'%+5.2f');
                        end
                        for i=1:obj.stem_par.p
                            for j=i+1:obj.stem_par.p
                                if not(isempty(obj.stem_EM_result.varcov))
                                    output{i+1,j+1}=[num2str(obj.stem_par.v_z(i,j),'%+05.2f'),' (Std ',num2str(sqrt(obj.stem_EM_result.varcov(counter,counter)),'%03.2f'),')'];
                                    output{j+1,i+1}=output{i+1,j+1};
                                else
                                    output{i+1,j+1}=num2str(obj.stem_par.v_z(i,j),'%+05.2f');
                                    output{j+1,i+1}=output{i+1,j+1};
                                end
                                counter=counter+1;
                            end
                        end
                        disp(output);
                    end
                end
            else
                disp('The model has not been estimated yet. Use the method print of the class stem_data to print data information.');
            end            
        end
        
        function set_system_size(obj)
            %DESCRIPTION: evaluate the minimum N after which it is faster to evaluate only the diagonal of a matrix product instead of the full matrix
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the system_size property is updated
            
            dim=100;
            t1=0;
            t2=1;
            while t1<t2
                m1=randn(dim);
                m2=randn(dim);
                tic
                diag(m1*m2);
                t1=toc;
                tic
                h2=zeros(dim,1);
                for i=1:size(m1)
                    h2(i)=m1(i,:)*m2(:,i);
                end
                t2=toc;
                dim=dim+50;
            end
            obj.system_size=dim;
        end
        
        function set_logL(obj)
            %DESCRIPTION: evaluate the observed data log-likelihood
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the logL property of the object stem_EM_result is updated
            disp('Log-Likelihood computation...');
            st_kalman=stem_kalman(obj);
            [st_kalmansmoother_result,~,~,~,~,~,~,~,~] = st_kalman.smoother(1);
            obj.stem_EM_result.logL=st_kalmansmoother_result.logL;
            disp('Log-Likelihood computation ended.');
        end
        
        function set_varcov(obj) 
            %DESCRIPTION: evaluate the variance-covariance matrix of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the varcov property of the object stem_EM_result is updated    
            
            %parameter order: beta,sigma_eps,alpha_bp,theta_b,v_b,alpha_p,theta_p,v_p,G,sigma_eta,[alpha_z,theta_z,v_z]
            
            if obj.estimated==0
                error('The model has not been estimated yet');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  dimension estimation   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            N=obj.N;
            Np=obj.Np;
            T=obj.T;
            dim=obj.dim;            
            
            data=obj.stem_data;
            par=obj.stem_par;
            q=par.q;
            p=par.p;
            k=par.k;

            n_beta=0;
            n_bp_alpha=0;
            n_bp_theta=0;
            n_bp_v=0;
            n_p_alpha=0;            
            n_p_theta=0;
            n_p_v=0;
            n_time_G=0;
            n_time_s2e=0;
            n_z_alpha=0;
            
            if not(isempty(data.X_bp))
                n_bp_alpha=2*q;
                if par.pixel_correlated
                    n_bp_v=q*(q-1)/2; %v_bp extra-diagonal and univoc
                    n_bp_theta=1;
                else
                    n_bp_theta=q;
                    n_bp_v=0;
                end
            end
            
            if not(isempty(data.X_beta))
                n_beta=length(par.beta);
            end
            
            if p>0
                if not(obj.stem_data.model_type==1)
                    if par.time_diagonal
                        n_time_G=p;
                        n_time_s2e=p;
                    else
                        n_time_G=p^2;
                        n_time_s2e=(p+1)*p/2;
                    end
                    n_time=n_time_G+n_time_s2e;
                    n_z_alpha=0;
                else
                    n_time_G=p;
                    n_z_alpha=p;
                    n_z_theta=1;
                    n_z_v=(p-1)*p/2;
                    n_time_s2e=0;
                    n_time=n_time_G+n_z_alpha+n_z_theta+n_z_v;
                end
            else
                n_time=0;
            end
            
            if not(isempty(data.X_p))
                n_p_alpha=q*k;
                n_p_v=q*(q-1)/2*k;
                n_p_theta=k;
            end
            n_eps=size(par.sigma_eps,1);
            
            n_psi=n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+n_p_v+n_time;
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  extraction of the useful data   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if p>0
                %kalman filter
                st_kalman=stem_kalman(obj);
                compute_logL=0;
                enable_varcov_computation=1;
                [st_kalmanfilter_result,sigma_eps,sigma_W_b,sigma_W_p,sigma_Z,sigma_eta,G_tilde_diag,sigma_geo,aj_bp,aj_p,aj_z,M] = st_kalman.filter(compute_logL,enable_varcov_computation);
                rr=size(sigma_Z,1);
                if not(data.model_type==1)
                    G=par.G;
                else
                    G=sparse(1:length(G_tilde_diag),1:length(G_tilde_diag),G_tilde_diag,length(G_tilde_diag),length(G_tilde_diag));
                end
            else
                [sigma_eps,sigma_W_b,sigma_W_p,sigma_geo,~,sigma_eta,~,aj_bp,aj_p,aj_z,M] = obj.get_sigma();
                st_kalmanfilter_result=stem_kalmanfilter_result([],[],[],[],[],[],[]);
                rr=0;
            end            
            %J=st_kalmanfilter_result.J(:,:,2:end); %J for t=0 is not considered
            J=st_kalmanfilter_result.J;
            zk_f=st_kalmanfilter_result.zk_f;
            Pk_f=st_kalmanfilter_result.Pk_f;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  derivatives allocation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Derivatives allocation...');
            d_e_Lt=cell(n_psi,1);
            if not(isempty(data.X_z))
                d_Z=cell(n_psi,1);
                d_J_Lt=cell(n_psi,1);
                d_P=cell(n_psi,1);
            end
            d_St_Lt=cell(n_psi,1);
            d_Sgeo=cell(n_psi,1);
            
            if not(isempty(data.X_beta))
                d_beta=cell(n_psi,1);
            end
            d_Seps=cell(n_psi,1);
            if not(isempty(data.X_bp))
                d_alpha_bp=cell(n_psi,n_bp_alpha);
                d_M_Sb=cell(n_psi,1);
            end
            if not(isempty(data.X_p))
                d_alpha_p=cell(n_psi,n_p_alpha,k);
                d_Sp=cell(n_psi,k);
            end
            if not(isempty(data.X_z))
                d_G=cell(n_psi,1);
                d_s2e=cell(n_psi,1);
                d_alpha_z=cell(n_psi,1);
                d_a_z_X=cell(n_psi,1);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  allocation of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:n_psi
                if not(isempty(data.X_beta))
                    d_beta{i}=sparse(n_beta,1);
                end
                d_Seps{i}=sparse(N,N);
                if not(isempty(data.X_bp))
                    for j=1:n_bp_alpha
                        d_alpha_bp{i,j}=sparse(N,1);
                    end
                    d_M_Sb{i}=sparse(N,N);
                end
                if not(isempty(data.X_p))
                    for j=1:k
                        d_Sp{i,j}=sparse(Np,Np);
                        for z=1:n_p_alpha
                            d_alpha_p{i,z,j}=sparse(N,1);
                        end
                    end
                end
                if not(isempty(data.X_z))
                    if obj.stem_data.model_type==1
                        for j=1:n_z_alpha
                            d_alpha_z{i,j}=sparse(N,1);
                        end
                    end
                    d_G{i}=sparse(rr,rr);
                    d_s2e{i}=sparse(rr,rr);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  initialization of time invariant derivatives  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Derivatives initialization...');
            
            blocks=[0 cumsum(dim)];
            
            %d_beta
            if not(isempty(data.X_beta))
                for i=1:n_beta
                    d_beta{i}(i,1)=1;
                end
            end

            %d_Seps
            for i=1:n_eps
                Id=blocks(i)+1:blocks(i+1);
                d_Seps{n_beta+i}=sparse(Id,Id,ones(length(Id),1),N,N);
            end
            
            if not(isempty(data.X_bp))
                %d_alpha_bp
                for i=1:n_bp_alpha
                    [~,j_b] = obj.get_jbp(i);
                    d_alpha_bp{n_beta+n_eps+i,i}=j_b;
                end
                
                %d_M_Sb with respect to theta_bp
                d=stem_misc.M_apply(obj.stem_data.DistMat_b,M,'b');
                result=stem_misc.M_apply(sigma_W_b,M,'b');
                if not(par.pixel_correlated)&&(q>1)
                    for j=1:q
                        Id=[];
                        Jd=[];
                        elements=[];
                        for i=1:2
                            idx=blocks(j+(i-1)*q)+1:blocks(j+(i-1)*q+1);
                            result1=d(idx,idx);
                            result2=result(idx,idx);
                            if strcmp(par.correlation_type,'exponential')
                                result3=result1/(par.theta_b(j)^2).*result2;
                            elseif strcmp(par.correleation_type,'matern32')
                                result3=-sqrt(3)*result1/(par.theta_b(j)^2).*result2./(1+sqrt(3).*result1./par.theta_b(j))+result2.*sqrt(3).*result1./(par.theta_b(j)^2);
                            else
                                result3=(-sqrt(5)*result1/(par.theta_b(j)^2)-10/3*result1.^2/par.theta_b(j)^3).*result2./(1+sqrt(5).*result1./par.theta_b(j)+5/3*result1.^2/par.theta_b(j)^2)+result2.*sqrt(5).*result1./(par.theta_b(j)^2);
                            end
                            
                            L=find(result3);
                            [idx_I,idx_J]=ind2sub(size(result3),L);
                            Id=cat(1,Id,idx_I+blocks(j+(i-1)*2));
                            Jd=cat(1,Jd,idx_J+blocks(j+(i-1)*2));
                            elements=cat(1,elements,result3(L));
                        end
                        zero_density=(1-length(Id)/(N^2))*100;
                        if obj.tapering||zero_density>60
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+j}=sparse(Id,Jd,elements,N,N);
                        else
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+j}=zeros(N);
                            temp=sub2ind([N,N],Id,Jd);
                            d_M_Sb{n_beta+n_eps+n_bp_alpha+j}(temp)=elements;
                        end
                    end
                else
                    d_M_Sb{n_beta+n_eps+n_bp_alpha+1}=result.*d/(par.theta_b^2);
                end
                
                %d_M_Sb with respect to v_bp
                if par.pixel_correlated
                    z=1;
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            for h=1:2
                                idx1=blocks(i+(h-1)*q)+1:blocks(i+(h-1)*q+1);
                                idx2=blocks(j+(h-1)*q)+1:blocks(j+(h-1)*q+1);
                                result1=result(idx1,idx2);
                                result2=result1/par.v_b(i,j);
                                L=find(result2);
                                [idx_I,idx_J]=ind2sub(size(result2),L);
                                Id=cat(1,Id,idx_I+blocks(i+(h-1)*q));
                                Jd=cat(1,Jd,idx_J+blocks(j+(h-1)*q));
                                elements=cat(1,elements,result2(L));
                                Id=cat(1,Id,idx_J+blocks(i+(h-1)*q));
                                Jd=cat(1,Jd,idx_I+blocks(j+(h-1)*q));
                                elements=cat(1,elements,result2(L));
                            end
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=sparse(Id,Jd,elements,N,N);
                            else
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}=zeros(N);
                                temp=sub2ind([N,N],Id,Jd);
                                d_M_Sb{n_beta+n_eps+n_bp_alpha+n_bp_theta+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                else
                    %nothing as the matrix V_b is not estimated
                end
            end

            if not(isempty(data.X_p))
                %d_alpha_p
                counter=1;
                for z=1:k
                    for i=1:q
                        [~,j_p] = obj.get_jp(i,z);
                        d_alpha_p{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+counter,i,z}=j_p;
                        counter=counter+1;
                    end
                end
                   
                %d_Sp with respect to theta_p
                for z=1:k
                    if strcmp(par.correlation_type,'exponential')
                        d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+z,z}=sigma_W_p{z}.*obj.stem_data.DistMat_p/(par.theta_p(z)^2);
                    elseif strcmp(par.correlation_type,'matern32')
                        d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+z,z}=-sqrt(3)*obj.stem_data.DistMat_p/(par.theta_p(z)^2).*exp(-sqrt(3)*obj.stem_data.DistMat_p/par.theta_p(z))+sigma_W_p{z}.*sqrt(3).*obj.stem_data.DistMat_p/(par.theta_p(z)^2);
                    else
                        d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+z,z}=(-sqrt(5)*obj.stem_data.DistMat_p/par.theta_p(z)^2-10/3.*obj.stem_data.DistMat_p.^2/par.theta_p(z)^3).*exp(-sqrt(5)*obj.stem_data.DistMat_p/par.theta_p(k))+sigma_W_p{z}.*sqrt(5).*obj.stem_data.DistMat_p/(par.theta_p(z)^2);
                    end
                    if stem_misc.zero_density(d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+k,k})>60
                        d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+z,z}=sparse(d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+z,z});
                    end
                end
                
                %d_Sp with respect to v_p
                counter=1;
                for z=1:k
                    for j=1:q
                        for i=j+1:q
                            Id=[];
                            Jd=[];
                            elements=[];
                            %since the block is extra-diagonal it has two separated D_apply!
                            result1=sigma_W_p{k}(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                            result2=result1/par.v_p(i,j,z);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=cat(1,Id,idx_I+blocks(j));
                            Jd=cat(1,Jd,idx_J+blocks(i));
                            elements=cat(1,elements,result2(L));
                            Id=cat(1,Id,idx_J+blocks(i));
                            Jd=cat(1,Jd,idx_I+blocks(j));
                            elements=cat(1,elements,result2(L));
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+counter,z}=sparse(Id,Jd,elements,Np,Np);
                            else
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+counter,z}=zeros(Np);
                                temp=sub2ind([Np,Np],Id,Jd);
                                d_Sp{n_beta+n_eps+n_bp_alpha+n_bp_theta+n_bp_v+n_p_alpha+n_p_theta+counter,z}(temp)=elements;
                            end
                            counter=counter+1;
                        end
                    end
                end                
            end
            
            if not(isempty(data.X_z))
                if not(obj.stem_data.model_type==1)
                    %d_G
                    if par.time_diagonal
                        for i=1:rr
                            d_G{n_psi-n_time+i}(i,i)=1;
                        end
                    else
                        for i=1:rr^2
                            [a,b]=ind2sub([rr rr],i);
                            d_G{n_psi-n_time+i}(a,b)=1;
                        end
                    end
                    %d_s2e
                    if par.time_diagonal
                        for i=1:p
                            d_s2e{n_psi-n_time_s2e+i}(i,i)=1;
                        end
                    else
                        j=1;
                        for i=1:(rr*(rr+1))/2
                            temp=zeros(rr);
                            l=1;
                            for h=1:rr
                                for z=h:rr
                                    if i==l
                                        temp(h,z)=1;
                                        temp(z,h)=1;
                                    end
                                    l=l+1;
                                end
                            end
                            d_s2e{n_psi-n_time_s2e+i}=sparse(temp);
                            j=j+1;
                        end
                    end
                else
                    %d_G
                    dim=obj.stem_data.dim;
                    blocks=[0 cumsum(dim(1:p))];
                    for i=1:p
                        j_z=zeros(rr,1);
                        j_z(blocks(i)+1:blocks(i+1))=1;
                        d_G{n_psi-n_time+i}=sparse(1:rr,1:rr,j_z,rr,rr);
                    end
                    
                    %d_alpha_z
                    for i=1:n_z_alpha
                        [~,j_z] = obj.get_jz(i);
                        d_alpha_z{n_psi-n_time+n_time_G+n_time_s2e+i,i}=j_z;
                    end
                    
                    %d_s2e
                    if obj.stem_data.model_subtype==0
                        if strcmp(par.correlation_type,'exponential')
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply(sigma_eta.*obj.stem_data.DistMat_p/(par.theta_z^2),aj_z,'b');
                        elseif strcmp(par.correlation_type,'matern32')
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply(-sqrt(3)*obj.stem_data.DistMat_p/(par.theta_z^2).*exp(-sqrt(3)*obj.stem_data.DistMat_p/par.theta_z)+sigma_eta.*sqrt(3).*obj.stem_data.DistMat_p/(par.theta_z^2),aj_z,'b');
                        else
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply((-sqrt(5)*obj.stem_data.DistMat_p/par.theta_z^2-10/3.*obj.stem_data.DistMat_p.^2/par.theta_z^3).*exp(-sqrt(5)*obj.stem_data.DistMat_p/par.theta_z)+sigma_eta.*sqrt(5).*obj.stem_data.DistMat_z/(par.theta_z^2),aj_z,'b');
                        end
                    else
                        if strcmp(par.correlation_type,'exponential')
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply(sigma_eta.*obj.stem_data.DistMat_z/(par.theta_z^2),aj_z,'b');
                        elseif strcmp(par.correlation_type,'matern32')
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply(-sqrt(3)*obj.stem_data.DistMat_z/(par.theta_z^2).*exp(-sqrt(3)*obj.stem_data.DistMat_z/par.theta_z)+sigma_eta.*sqrt(3).*obj.stem_data.DistMat_z/(par.theta_z^2),aj_z,'b');
                        else
                            d_s2e{n_psi-n_time+p+p+1}=stem_misc.D_apply((-sqrt(5)*obj.stem_data.DistMat_z/par.theta_z^2-10/3.*obj.stem_data.DistMat_z.^2/par.theta_z^3).*exp(-sqrt(5)*obj.stem_data.DistMat_z/par.theta_z)+sigma_eta.*sqrt(5).*obj.stem_data.DistMat_z/(par.theta_z^2),aj_z,'b');
                        end
                    end
                    
                    %v_z
                    z=1;
                    for j=1:p
                        for i=j+1:p
                            Id=[];
                            Jd=[];
                            elements=[];
                            %since the block is extra-diagonal it has two separated D_apply!
                            result1=sigma_eta(blocks(j)+1:blocks(j+1),blocks(i)+1:blocks(i+1));
                            result2=stem_misc.D_apply(result1,aj_z(blocks(j)+1:blocks(j+1)),'l');
                            result2=stem_misc.D_apply(result2,aj_z(blocks(i)+1:blocks(i+1)),'r');
                            result2=result2/par.v_z(i,j);
                            L=find(result2);
                            [idx_I,idx_J]=ind2sub(size(result2),L);
                            Id=cat(1,Id,idx_I+blocks(j));
                            Jd=cat(1,Jd,idx_J+blocks(i));
                            elements=cat(1,elements,result2(L));
                            Id=cat(1,Id,idx_J+blocks(i));
                            Jd=cat(1,Jd,idx_I+blocks(j));
                            elements=cat(1,elements,result2(L));
                            zero_density=(1-length(Id)/(N^2))*100;
                            if obj.tapering||zero_density>60
                                d_s2e{n_psi-n_time+p+p+1+z}=sparse(Id,Jd,elements,Np,Np);
                            else
                                d_s2e{n_psi-n_time+p+p+1+z}=zeros(Np);
                                temp=sub2ind([Np,Np],Id,Jd);
                                d_s2e{n_psi-n_time+p+p+1+z}(temp)=elements;
                            end
                            z=z+1;
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  preliminary computations for time-invariant case  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if not(data.X_tv) 
                for i=1:n_psi
                    d_Sgeo{i}=sparse(N,N);
                end
                if not(isempty(data.X_bp))
                    a_bp_X=stem_misc.D_apply(data.X_bp(:,:,1),aj_bp,'l');
                    for i=1:n_psi
                        temp=zeros(N,1);
                        for j=1:n_bp_alpha
                            temp=temp+stem_misc.D_apply(data.X_bp(:,:,1),d_alpha_bp{i,j},'l');
                        end
                        d_a_bp_X=sparse(temp);
                        temp=stem_misc.M_apply(sigma_W_b,M,'b');
                        temp=stem_misc.D_apply(temp,d_a_bp_X,'l');
                        temp=stem_misc.D_apply(temp,a_bp_X,'r');
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                        
                        temp=stem_misc.D_apply(d_M_Sb{i},a_bp_X,'b');
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                        
                        temp=stem_misc.M_apply(sigma_W_b,M,'b');
                        temp=stem_misc.D_apply(temp,d_a_bp_X,'r');
                        temp=stem_misc.D_apply(temp,a_bp_X,'l');                        
                        d_Sgeo{i}=d_Sgeo{i}+temp;
                    end
                end
                if not(isempty(data.X_p))
                    for z=1:k
                        a_p_X=stem_misc.D_apply(data.X_p(:,:,1,z),aj_p(:,z),'l');
                        for i=1:n_psi
                            temp=zeros(N,1);
                            for j=1:q
                                temp=temp+stem_misc.D_apply(data.X_p(:,:,1,z),d_alpha_p{i,j,z},'l');
                            end
                            d_a_p_X=sparse(temp);
                            
                            temp=stem_misc.D_apply(sigma_W_p{z},d_a_p_X,'l');
                            temp=stem_misc.D_apply(temp,a_p_X,'r');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                            
                            temp=stem_misc.D_apply(d_Sp{i,z},a_p_X,'b');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                            
                            temp=stem_misc.D_apply(sigma_W_p{z},d_a_p_X,'r');
                            temp=stem_misc.D_apply(temp,a_p_X,'l');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                        end
                    end
                end
                for i=1:n_psi
                    d_Sgeo{i}=d_Sgeo{i}+d_Seps{i};
                end
                
                if not(isempty(data.X_z))
                    if (obj.stem_data.model_type==1)&&(obj.stem_data.model_subtype==0)
                        temp=data.X_z(:,:,1);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[data.X_z(:,:,1);zeros(N-size(data.X_z(:,:,1),1),size(data.X_z(:,:,1),2))];
                    end
                    a_z_X=stem_misc.D_apply(X_z_orlated,aj_z,'l'); 
                end
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Information matrix evaluation  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp('Hessian evaluation...');
            c0 = 0;
            IM=zeros(n_psi);
            tot=n_psi*(n_psi+1)/2*T;
            counter=1;
            
            for t=2:T
                if data.X_bp_tv
                    tBP=t-1;
                else
                    tBP=1;
                end
                if data.X_z_tv
                    tT=t-1;
                else
                    tT=1;
                end
                if data.X_beta_tv
                    tbeta=t-1;
                else
                    tbeta=1;
                end      
                if data.X_p_tv
                    tP=t-1;
                else
                    tP=1;
                end                  
                
                if data.X_tv
                    %compute sigma_geo in the time-variant case
                    if not(isempty(data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                    end
                    if not(isempty(data.X_p))
                        if isempty(data.X_bp)
                            if (obj.tapering)
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for z=1:k
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{z},data.X_p(:,1,tP,z),'b'),aj_p(:,z),'b');
                        end
                    end
                    if isempty(data.X_p)&&isempty(data.X_bp)
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    if not(isempty(data.X_bp))
                        sigma_geo=stem_misc.D_apply(stem_misc.D_apply(stem_misc.M_apply(sigma_W_b,M,'b'),data.X_bp(:,1,tBP),'b'),aj_bp,'b');
                    end
                    if not(isempty(data.X_p))
                        if not(isempty(sigma_geo))
                            if obj.tapering
                                sigma_geo=spalloc(size(sigma_W_p{1},1),size(sigma_W_p{1},1),nnz(sigma_W_p{1}));
                            else
                                sigma_geo=zeros(N);
                            end
                        end
                        for z=1:k
                            sigma_geo=sigma_geo+stem_misc.D_apply(stem_misc.D_apply(sigma_W_p{z},data.X_p(:,1,tP,z),'b'),aj_p(:,z),'b');
                        end
                    end
                    if not(isempty(sigma_geo))
                        sigma_geo=sigma_eps;
                    else
                        sigma_geo=sigma_geo+sigma_eps;
                    end
                    
                    %compute d_Sgeo in the time-variant case
                    for i=1:n_psi
                        d_Sgeo{i}=sparse(N,N);
                    end
                    if not(isempty(data.X_bp))
                        a_bp_X=stem_misc.D_apply(data.X_bp(:,:,tBP),aj_bp,'l');
                        for i=1:n_psi
                            temp=zeros(N,1);
                            for j=1:n_bp_alpha
                                temp=temp+stem_misc.D_apply(data.X_bp(:,:,tBP),d_alpha_bp{i,j},'l');
                            end
                            d_a_bp_X=sparse(temp);
                            temp=stem_misc.M_apply(sigma_W_b,M,'b');
                            temp=stem_misc.D_apply(temp,d_a_bp_X,'l');
                            temp=stem_misc.D_apply(temp,a_bp_X,'r');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                            
                            temp=stem_misc.D_apply(d_M_Sb{i},a_bp_X,'b');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                            
                            temp=stem_misc.M_apply(sigma_W_b,M,'b');
                            temp=stem_misc.D_apply(temp,d_a_bp_X,'r');
                            temp=stem_misc.D_apply(temp,a_bp_X,'l');
                            d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                        end
                    end
                    if not(isempty(data.X_p))
                        for z=1:k
                            a_p_X=stem_misc.D_apply(data.X_p(:,:,tP,z),aj_p(:,z),'l');
                            for i=1:n_psi
                                temp=zeros(N,1);
                                for j=1:q
                                    temp=temp+stem_misc.D_apply(data.X_p(:,:,tP,z),d_alpha_p{i,j,z},'l');
                                end
                                d_a_p_X=sparse(temp);
                                
                                temp=stem_misc.D_apply(sigma_W_p{z},d_a_p_X,'l');
                                temp=stem_misc.D_apply(temp,a_p_X,'r');
                                d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                                
                                temp=stem_misc.D_apply(d_Sp{i,z},a_p_X,'b');
                                d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                                
                                temp=stem_misc.D_apply(sigma_W_p{z},d_a_p_X,'r');
                                temp=stem_misc.D_apply(temp,a_p_X,'l');
                                d_Sgeo{i}=d_Sgeo{i}+stem_misc.sparseif(temp,60);
                            end
                        end
                    end
                    for i=1:n_psi
                        d_Sgeo{i}=d_Sgeo{i}+d_Seps{i};
                    end
                    
                    if not(isempty(data.X_z))
                        if (obj.stem_data.model_type==1)&&(obj.stem_data.model_subtype==0)
                            temp=data.X_z(:,:,tT);
                            temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                            X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                        else
                            X_z_orlated=[data.X_z(:,:,1);zeros(N-size(data.X_z(:,:,tT),1),size(data.X_z(:,:,tT),2))];
                        end
                        a_z_X=stem_misc.D_apply(X_z_orlated,aj_z,'l');
                    end
                end
                
                Lt=not(isnan(data.Y(:,t-1)));
                if t>2
                    Lt_lag=not(isnan(data.Y(:,t-2)));
                end
                if not(isempty(data.X_z))
                    if (obj.stem_data.model_type==1)&&(obj.stem_data.model_subtype==0)
                        temp=data.X_z(:,:,tT);
                        temp=sparse(1:length(temp),1:length(temp),temp,length(temp),length(temp));
                        X_z_orlated=cat(1,temp,zeros(N-size(temp,1),size(temp,2)));
                    else
                        X_z_orlated=[data.X_z(:,:,tT);zeros(N-size(data.X_z(:,:,tT),1),size(data.X_z(:,:,tT),2))];
                    end
                    X_z_orlated=stem_misc.D_apply(X_z_orlated,aj_z,'l');                    
                    
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                        e_t_Lt=data.Y(Lt,t-1)-X_beta_orlated(Lt,:)*par.beta-X_z_orlated(Lt,:)*zk_f(:,t);
                    else
                        e_t_Lt=data.Y(Lt,t-1)-X_z_orlated(Lt,:)*zk_f(:,t);
                    end                    
                    
                    sigma_t_Lt=X_z_orlated(Lt,:)*Pk_f(:,:,t)*X_z_orlated(Lt,:)'+sigma_geo(Lt,Lt);
                    
                    d_Sz=cell(n_psi,1);
                    for i=1:n_psi
                        temp=zeros(N,size(X_z_orlated,2));
                        for j=1:n_z_alpha
                            temp=temp+stem_misc.D_apply(X_z_orlated,d_alpha_z{i,j},'l');
                        end
                        d_a_z_X{i}=sparse(temp);
                        
                        d_Sz{i}=sparse(N,N);
                    end
                else
                    if n_beta>0
                        X_beta_orlated=data.X_beta(:,:,tbeta);
                        X_beta_orlated=cat(1,X_beta_orlated,zeros(N-size(X_beta_orlated,1),size(X_beta_orlated,2)));
                        e_t_Lt=data.Y(Lt,t-1)-X_beta_orlated(Lt,:)*par.beta;
                    else
                        e_t_Lt=data.Y(Lt,t-1);
                    end
                    sigma_t_Lt=sigma_geo(Lt,Lt);
                end
                
                sigma_t_Lt=stem_misc.sparseif(sigma_t_Lt,60);
                
                for i=1:n_psi
                    if t==2
                        if not(isempty(data.X_z))
                            d_P{i}=d_s2e{i};
                            
                            temp=d_a_z_X{i}*Pk_f(:,:,t)*a_z_X';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            temp=a_z_X*d_P{i}*a_z_X';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            temp=a_z_X*Pk_f(:,:,t)*d_a_z_X{i}';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            
                            d_St_Lt{i}=d_Sz{i}(Lt,Lt)+d_Sgeo{i}(Lt,Lt);
                            %verificare sistema lineare nel caso di matrici sparse!
                            d_J_Lt{i}=(d_G{i}*Pk_f(:,:,t)*a_z_X(Lt,:)'+G*d_P{i}*a_z_X(Lt,:)'+G*Pk_f(:,:,t)*d_a_z_X{i}(Lt,:)'-J(:,Lt,t)*d_St_Lt{i})/sigma_t_Lt;
                            d_Z{i}=d_G{i}*zk_f(:,t-1);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        end
                    else
                        if not(isempty(data.X_z))
                            d_P{i}=d_G{i}*Pk_f(:,:,t-1)*G+G*d_P_lag{i}*G'+G*Pk_f(:,:,t-1)*d_G{i}'+d_s2e{i}-d_J_Lt_lag{i}*sigma_t_Lt_lag*J(:,Lt_lag,t-1)'-J(:,Lt_lag,t-1)*d_St_Lt_lag{i}*J(:,Lt_lag,t-1)'-J(:,Lt_lag,t-1)*sigma_t_Lt_lag*d_J_Lt_lag{i}';
                            
                            temp=d_a_z_X{i}*Pk_f(:,:,t)*a_z_X';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            temp=a_z_X*d_P{i}*a_z_X';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            temp=a_z_X*Pk_f(:,:,t)*d_a_z_X{i}';
                            d_Sz{i}=d_Sz{i}+stem_misc.sparseif(temp,60);
                            
                            d_St_Lt{i}=d_Sz{i}(Lt,Lt)+d_Sgeo{i}(Lt,Lt);
                            %verificare sistema lineare nel caso di matrici sparse!
                            d_J_Lt{i}=(d_G{i}*Pk_f(:,:,t)*a_z_X(Lt,:)'+G*d_P{i}*a_z_X(Lt,:)'+G*Pk_f(:,:,t)*d_a_z_X{i}(Lt,:)'-J(:,Lt,t)*d_St_Lt{i})/sigma_t_Lt;
                            d_Z{i}=d_G{i}*zk_f(:,t-1)+G*d_Z_lag{i}+d_J_Lt_lag{i}*e_t_Lt_lag+J(:,Lt_lag,t-1)*d_e_Lt_lag{i};
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i}-d_a_z_X{i}(Lt,:)*zk_f(:,t)-a_z_X(Lt,:)*d_Z{i};
                            else
                                d_e_Lt{i}=-d_a_z_X{i}(Lt,:)*zk_f(:,t)-a_z_X(Lt,:)*d_Z{i};
                            end
                        else
                            d_St_Lt{i}=d_Sgeo{i}(Lt,Lt);
                            if n_beta>0
                                d_e_Lt{i}=-X_beta_orlated(Lt,:)*d_beta{i};
                            else
                                d_e_Lt{i}=zeros(sum(Lt),1);
                            end
                        end
                    end
                end
                
                temp0=cell(n_psi,1);
                temp1=cell(n_psi,1);
                
                for i=1:n_psi
                    if issparse(sigma_t_Lt)
                        r = symamd(sigma_t_Lt);
                        c=chol(sigma_t_Lt(r,r));
                        temp0{i}(r,:)=stem_misc.chol_solve(c,d_e_Lt{i}(r,:));
                    else
                        c=chol(sigma_t_Lt);
                        temp0{i}=stem_misc.chol_solve(c,d_e_Lt{i});
                    end
                    
                    d_St_i_Lt=d_St_Lt{i};
                    if nnz(d_St_i_Lt)>0
                        if issparse(sigma_t_Lt)
                            temp1{i}(r,:)=full(stem_misc.chol_solve(c,d_St_i_Lt(r,:)));
                        else
                            temp1{i}=full(stem_misc.chol_solve(c,d_St_i_Lt));
                        end
                    else
                        temp1{i}=spalloc(size(sigma_t_Lt,1),size(sigma_t_Lt,2),0);
                    end
                end
                 
                for i=1:n_psi
                    blocks=0:50:size(temp1{i},1);
                    if not(blocks(end)==size(temp1{i},1))
                        blocks=cat(2,blocks,size(temp1{i},1));
                    end
                    nnz_temp1_i=nnz(temp1{i});
                    
                    for j=i:n_psi
                        IM(i,j)=IM(i,j)+d_e_Lt{i}'*temp0{j};
                        if (nnz_temp1_i>0)&&(nnz(temp1{j})>0)
                            sumtrace=0;
                            for z=1:length(blocks)-1
                                idx=blocks(z)+1:blocks(z+1);
                                sumtrace=sumtrace+trace(temp1{i}(idx,:)*temp1{j}(:,idx));
                            end
                            IM(i,j)=IM(i,j)+0.5*sumtrace+0.25*trace(temp1{i})*trace(temp1{j});
                        end
                        counter=counter+1;
                        if (mod(counter,100)==0)||(counter==tot)
                            if mod(round(counter/tot*100),1)==0&&c0~=round(counter/tot*100)
                                c0 = round(counter/tot*100);
                                disp(['Hessian evaluation: ',num2str(round(counter/tot*100)),'% completed']);
                            end
                        end
                    end
                end

                if not(isempty(data.X_z))
                    d_P_lag=d_P;
                    d_J_Lt_lag=d_J_Lt;
                    d_Z_lag=d_Z;
                end

                d_e_Lt_lag=d_e_Lt;
                e_t_Lt_lag=e_t_Lt;
                d_St_Lt_lag=d_St_Lt;
                sigma_t_Lt_lag=sigma_t_Lt;
                if data.X_tv
                    sigma_geo=[];
                end
            end
            IM=IM+triu(IM,1)';
            obj.stem_EM_result.varcov=inv(IM);
        end      
        
        function set_initial_values(obj,stem_par)
            %DESCRIPTION: set the initial values of the model parameters
            %
            %INPUT
            %
            %obj - [stem_model object] (1x1)
            %
            %OUTPUT
            %
            %none: the stem_par_initial property is updated                
            if not(isa(stem_par,'stem_par'))
                error('The input argument must be of class stem_par');
            end
            obj.stem_par_initial=stem_par;
            if obj.stem_par.time_diagonal&&not(obj.stem_data.model_type==1)
                obj.stem_par_initial.G=diag(diag(obj.stem_par_initial.G));
                obj.stem_par_initial.sigma_eta=diag(diag(obj.stem_par_initial.sigma_eta));
            end
        end
         
        %Export functions. Useful to avoid access to the properties of nested objects
        function N = N(obj)
            N=obj.stem_data.N();
        end
        
        function Nb = Nb(obj)
            Nb=obj.stem_data.Nb();
        end
        
        function Np = Np(obj)
            Np=obj.stem_data.Np();
        end
        
        function T = T(obj)
            T=obj.stem_data.T();
        end
        
        function nvar=nvar(obj)
            nvar=obj.stem_data.nvar();
        end
        
        function dim=dim(obj)
            dim=obj.stem_data.dim();
        end
        
        %Initial values estimation functions (only beta at the moment)
        function [beta0] = get_beta0(obj)
            if obj.stem_par.n_beta>0
                N = size(obj.stem_data.X_beta,1);
                y = obj.stem_data.Y(1:N,:);
                y=y(:);
                T = obj.T;
                x = zeros(N*T,size(obj.stem_data.X_beta,2));
                
                for t=1:T
                    if size(obj.stem_data.X_beta,3)==T
                        tT=t;
                    else
                        tT=1;
                    end
                    x((t-1)*N+1:t*N,:)=obj.stem_data.X_beta(:,:,tT);
                end
                
                L=not(isnan(y));
                beta0 = (x(L,:)'*x(L,:))\x(L,:)'*y(L);
            else
                disp('WARNING: the model does not include data to estimate beta');
                beta0=[];
            end
        end
      
        %Class set functions 
        function set.stem_data(obj,stem_data)
            if isa(stem_data,'stem_data')
                obj.stem_data=stem_data;
            else
                error('The argument must be of class stem_data');
            end
        end
        
        function set.stem_par(obj,stem_par)
            if not(isa(stem_par,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            end
            obj.stem_par=stem_par;
        end        
       
        function set.stem_par_initial(obj,stem_par_initial)
            if not(isa(stem_par_initial,'stem_par'))
                error('stem_par_initial must be of class stem_par');
            end
            obj.stem_par_initial=stem_par_initial;
        end
        
    end
    
end
