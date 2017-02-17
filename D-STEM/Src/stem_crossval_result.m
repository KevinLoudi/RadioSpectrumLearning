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



classdef stem_crossval_result < handle
    
    %CONSTANTS
    %
    %dN - the number of cross-validation sites
    %T  - the number of time steps
    
    properties
        stem_krig_result=[];        %[strm_krig_result object] (1x1)  the prediction over the cross-validation sites is obtained through kriging
        
        y_back=[];                  %[double]                  (dNxT) the original data (back-transformed) of the cross-validation sites
        y_hat_back=[];              %[double]                  (dNxT) the estimated data (back-transformed) for the cross-validation sites
        res=[];                     %[double]                  (dNxT) cross-validation residuals
        res_back=[];                %[double]                  (dNxT) back-transformed cross-validation residuals (if the original data are log-transformed and/or standardized)
        
        mse=[];                     %[double >=0]              (dNx1) mean squared error for each cross-validation site
        mse_time=[];                %[double >=0]              (dNx1) mean squared error for each time step
        relative_mse=[];            %[double >=0]              (dNx1) relative mean squared error for each cross-validation site
        relative_mse_time=[];       %[double >=0]              (dNx1) relative mean squared error for each time step
        min_distance=[];            %[double >=0]              (dNx1) the distance from each cross-validation site to the nearest non cross-validation site
    end
    
    methods
        function obj = stem_crossval_result()
            %DESCRIPTION: object constructor
            %
            %INPUT
            %no input required
            %
            %
            %OUTPUT
            %obj           - [stem_crossval_result object] (1x1)
        end
        
        function set.stem_krig_result(obj,stem_krig_result)
            if not(isa(stem_krig_result,'stem_krig_result'))
                error('stem_krig_result must be of class stem_krig_result');
            end
            obj.stem_krig_result=stem_krig_result;
        end
    end
end