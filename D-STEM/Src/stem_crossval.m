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



classdef stem_crossval < handle
    
    %CONSTANTS
    %
    %dq - number of cross-validation variables
    
    properties
        variable_name={};           %[string]                      {dqx1} the names of the cross-validation variables
        type={};                    %[string]                      {dqx1} 'point': cross-validation is on point-data
        indices={};                 %[integer >0]                  {dq}x(dNx1) the indices of the cross-validation sites 
        
        stem_varset={};             %[stem_varset objects]         {dqx1} the subset of data used for cross-validation
        stem_gridlist={};           %[stem_gridlist objects]       {dqx1} this object is needed for storing the coordinates of the cross-validation sites
        stem_crossval_result={};    %[stem_crossval_result object] {dqx1} the objects including the cross-validation results for each variable
    end

    methods
        function obj = stem_crossval(variable_name,type,indices)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %variable_name - [string]               {dqx1} the names of the cross-validation variables
            %type          - [string]               {dqx1} 'point': cross-validation is on point-data
            %indices       - [integer >0]           {dq}x(dNx1) the indices of the cross-validation sites for each variable
            %
            %OUTPUT
            %obj           - [stem_crossval object] (1x1)    
            
            if nargin<3
                error('Not enough input arguments');
            end
            if not(iscell(variable_name))
                error('variable_name must be a cell array');
            end
            if not(iscell(type))
                error('type must be a cell array');
            end
            if not(iscell(indices))
                error('indices must be a cell array');
            end
            if not(length(variable_name)==length(type) && length(type)==length(indices))
                error('variable name, type and indices must be cell array of the same length');
            end
            obj.variable_name=variable_name;
            obj.type=type;            
            obj.indices=indices;
        end
       
        %Class set methods
        function set.type(obj,type)
            for i=1:length(type)
                if not(strcmp(type{i},'point'))
                    error('The cross-validation is available only for point variables. The type must be equal to ''point.''');
                end
            end
            obj.type=type;
        end
        
        function set.indices(obj,indices)
            for i=1:length(indices)
                if min(indices{i})<1
                    error('The indices vector cannot contain negative values');
                end
            end
            obj.indices=indices;
        end            
    end
    
end