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

classdef stem_kalmanfilter_result < handle
    properties
        zk_f   = [];    %[double]     (pxT+1)    the filtered state
        zk_u   = [];    %[double]     (pxT+1)    the updated state
        Pk_f   = [];    %[double]     (pxpxT+1)  variance-covariance matrix of the filtered state
        Pk_u   = [];    %[double]     (pxpxT+1)  variance-covariance matrix of the updated state
        J_last = [];    %[double]     (pxN)      innovation vector at time t=T
        J      = [];    %[double]     (pxNxT+1)  innovation vector from time t=0 to time t=T
        logL   = [];    %[double]     (1x1)      observed-data log-likelihood
    end
    
    methods
        function obj = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL)
            %DESCRIPTION: constructor of the class stem_kalmanfilter_result
            %
            %INPUT 
            %See the class properties
            %
            %OUTPUT
            %obj             - [stem_kalmanfilter_result object]   (1x1) stem_kalmanfilter_result object  
            
            obj.zk_f=zk_f;
            obj.zk_u=zk_u;
            obj.Pk_f=Pk_f;
            obj.Pk_u=Pk_u;
            obj.J_last=J_last;
            obj.J=J;
            obj.logL=logL;
        end
    end
end