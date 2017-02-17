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



classdef stem_datestamp < handle
    
    properties
        date_start=[];      %[integer>0]  (1x1) the date of the first time step. It can be a date in the Matlab format (the output of the datenum function) or a numeric index
        date_end=[];        %[integer>0]  (1x1) the date of the last  time step. It can be a date in the Matlab format (the output of the datenum function) or a numeric index
        T=[];               %[integer>0]  (1x1) the total number of time steps
    end
    
    properties (SetAccess = private)
        stamp=[];           %[integer >0] (Tx1) all the date stamps
    end
    
    methods
        function obj = stem_datestamp(date_start,date_end,T)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %date_start    - [string|integer>0]     (1x1) the date related to the first time step. It can be a string in the format dd-mm-yyyy HH:MM or an integer index 
            %date_end      - [string|integer>0]     (1x1) the date related to the last  time step. It can be a string in the format dd-mm-yyyy HH:MM or an integer index
            %T             - [integer>0]            (1x1) the total number of time steps
            %
            %OUTPUT
            %obj           - [stem_datestamp object](1x1)  
            
            if nargin<3
                error('Not enough input arguments');
            end
            if isnumeric(date_start)
                obj.date_start=date_start;
            else
                obj.date_start=datenum(date_start,'dd-mm-yyyy HH:MM');
            end
            if isnumeric(date_end)
                obj.date_end=date_end;
            else
                obj.date_end=datenum(date_end,'dd-mm-yyyy HH:MM');
            end
            if obj.date_end<obj.date_start
                error('date_start cannot be higher than date_end');
            end
            obj.T=T;
            obj.stamp=obj.date_start:(obj.date_end-obj.date_start+1)/obj.T:obj.date_end;
        end
        
        function subset_stamps(obj,indices)
            %DESCRIPTION: removes a subset of the date stamps
            %
            %INPUT
            %obj           - [stem_datestamp object]    (1x1)  stem_datestamp object
            %indices       - [integer>0]                (dTx1) the indices of the temporal steps to keep
            %
            %OUTPUT
            %none: the properties of the object are updated
            
            obj.stamp=obj.stamp(indices);
            obj.date_start=min(obj.stamp);
            obj.date_end=max(obj.stamp);
            obj.T=length(obj.stamp);
        end
        
        function average_stamps(obj,indices)
            %DESCRIPTION: averages the date stamps. This method is used when the time_average method of the class stem_data is called
            %
            %INPUT
            %obj           - [stem_datestamp object]    (1x1) stem_datestamp object
            %indices       - [integer>0]                (1x1) the indices of the reference temporal steps. The average is evaluated on date stamps between indices(i)+1 and indices(i+1) included
            %
            %OUTPUT
            %none: the properties of the object are updated   
            
            stamp_temp=[];
            for i=1:length(indices)-1
                stamp_temp(i)=mean(obj.stamp(indices(i)+1):obj.stamp(indices(i+1)));
            end
            obj.stamp=stamp_temp;
            obj.date_start=min(obj.stamp);
            obj.date_end=max(obj.stamp);
            obj.T=length(obj.stamp);
        end
    end
end