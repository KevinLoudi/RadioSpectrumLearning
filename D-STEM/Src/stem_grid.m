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

classdef stem_grid
    
    %CONSTANTS
    %
    %ni_g - the number of point sites for the i-th variable
    %ni_r - the number of pixel sites for the i-th variable
    
    properties
        coordinate=[];              %[double]    (ni_g|ni_rx2)  matrix of spatial coordinates composed of either x and y vectors or latitude and longitude vectors
        unit='';                    %[string]    (1x1)          unit of measure of the coordinates. 'deg': degree; 'km': kilometers; 'm': meters
        grid_type=[];               %[string]    (1x1)          'sparse': sparse spatial locations; 'regular': spatial locations at fixed intervals in space
        grid_size=[];               %[integer >0](2x1)          if grid_type='regular', then grid_size is the number of rows and columns of the grid
        site_type=[];               %[string]    (1x1)          'point': the spatial coordinates relate to point data; 'pixel': the spatial coordinates related to pixel data
        pixel_shape='';             %[string]    (1x1)          'square': square pixels; 'rectangular': rectangular pixels
        pixel_side_w=[];            %[double]    (1x1)          the width of the pixel (in the same unit of measure of the unit property)
        pixel_side_h=[];            %[double]    (1x1)          the height of the pixel (in the same unit of measure of the unit property)
        duplicated_sites=[];        %[integer]   (hx1)          indices of the duplicated spatial locations in the coordinate property 
    end
    
    properties (Dependent, SetAccess = private)
       box=[];                      %[double]    (4x1) the bounding box of the geographic area covered by the grid [lat_min,lat_max,lon_min,lon_max] 
    end
    
    methods
        function obj = stem_grid(coordinate,unit,grid_type,site_type,grid_size,pixel_shape,pixel_side_w,pixel_side_h)
            %DESCRIPTION: is the constructor of the class stem_data
            %
            %INPUT
            %
            %coordinate         - [double]              (ni_g|ni_rx2)  matrix of spatial coordinates composed of either x and y vectors or latitude and longitude vectors
            %unit='';           - [string]              (1x1)          unit of measure of the coordinates. 'deg': degree; 'km': kilometers; 'm': meters
            %grid_type=[];      - [string]              (1x1)          'sparse': sparse spatial locations; 'regular': spatial locations at fixed intervals in space
            %grid_size=[];      - [integer >0]          (2x1)          if grid_type='regular', then grid_size is the number of rows and columns of the grid
            %site_type=[];      - [string]              (1x1)          'point': the spatial coordinates relate to point data; 'pixel': the spatial coordinates related to pixel data
            %pixel_shape='';    - [string]              (1x1)          'square': square pixels; 'rectangular': rectangular pixels
            %pixel_side_w=[];   - [double]              (1x1)          the width of the pixel (in the same unit of measure of the unit property)
            %pixel_side_h=[];   - [double]              (1x1)          the height of the pixel (in the same unit of measure of the unit property)
            %
            %OUTPUT
            %obj                - [stem_grid object]    (1x1)
            if nargin<4
                error('Not enough input arguments');
            end
            if strcmp(site_type,'pixel')&&(nargin<8)
                error('pixel_shape, pixel_side_w and pixel_side_h must be provided');
            end
            if nargin==4
                grid_size=[];
            end
            obj.grid_type=grid_type; %the other of this two lines is important and cannot change
            obj.coordinate=coordinate;
            obj.unit=unit;
            obj.site_type=site_type;
            obj.grid_size=grid_size;
            if nargin>5
                if not(strcmp(site_type,'pixel'))
                    warning('pixel_shape in ignored');
                else
                    obj.pixel_shape=pixel_shape;
                end
            end
            if nargin==7
                error('Also pixel_side_h must be provided');
            end
            if nargin>=8
                if not(strcmp(site_type,'pixel'))
                    warning('pixel_side is ignored');
                else
                    obj.pixel_side_w=pixel_side_w;
                    obj.pixel_side_h=pixel_side_h;
                end                
            end
        end
        
        % Class set methods
        function box = get.box(obj)
            if not(isempty(obj.coordinate))
                box(1)=min(obj.coordinate(:,1));
                box(2)=max(obj.coordinate(:,1));
                box(3)=min(obj.coordinate(:,2));
                box(4)=max(obj.coordinate(:,2));
            else
                box=0;
            end
        end
        
        function obj = set.coordinate(obj,coordinate)
            if not(size(coordinate,2)==2)
                error('coordinate must be a Nx2 matrix');
            end
            
            if not(strcmp(obj.grid_type,'regular'))
                if length(coordinate)<5000
                    obj.duplicated_sites=[];
                    for i=1:length(coordinate)-1
                        temp=coordinate(i,:);
                        temp2=coordinate((i+1):end,:);
                        temp_lat=temp2(:,1);
                        temp_lon=temp2(:,2);
                        a=temp_lat==temp(1);
                        b=temp_lon==temp(2);
                        c=a&b;
                        if sum(c)>0
                            obj.duplicated_sites=[obj.duplicated_sites;find(c,1)+i];
                            disp(['WARNING: coordinate ',num2str(i),' equal to coordinate ',num2str(find(c,1)+i)]);
                        end
                    end
                else
                    disp(['WARNING: too many spatial locations. The test for duplicate spatial locations is skipped']);
                end
            end

            obj.coordinate=coordinate;
        end
        
        function obj = set.unit(obj,unit)
            if not(strcmp(unit,'deg') || strcmp(unit,'m') || strcmp(unit,'km'))
                error('unit must be ''deg'' or ''m'' or ''km''');
            end
            obj.unit=unit;
        end
        
        function obj = set.grid_type(obj,grid_type)
            if not(strcmp(grid_type,'regular')||strcmp(grid_type,'sparse'))
                error('The grid type must be either regular or sparse');
            end
            obj.grid_type=grid_type;
        end
        
        function obj = set.site_type(obj,site_type)
            if not(strcmp(site_type,'point')||strcmp(site_type,'pixel'))
                error('The grid type must be either point or pixel');
            end
            obj.site_type=site_type;            
        end
        
        function obj = set.grid_size(obj,grid_size)
            if (not(isempty(grid_size)))&&(strcmp(obj.grid_type,'sparse'))
                error('The grid size must be provided only for regular grid');
            end
            if not(isempty(grid_size))
                if not(isvector(grid_size))
                    error('The grid_size must be a 2x1 vector');
                end
                if not(grid_size(1)*grid_size(2)==size(obj.coordinate,1))
                    error('The grid size is not compatible with the grid');
                end
            end
            obj.grid_size=grid_size;
        end
        
        function obj = set.pixel_shape(obj,pixel_shape)
            if not(strcmp(pixel_shape,'square') || strcmp(pixel_shape,'rectangular'))
                error('pixel_shape can be only ''square'' or ''rectangular''');
            end
            obj.pixel_shape=pixel_shape;
        end
        
        function obj = set.pixel_side_w(obj,pixel_side_w)
            if pixel_side_w<0
                error('pixel_side_w must be >0');
            end
            obj.pixel_side_w=pixel_side_w;
        end
        
        function obj = set.pixel_side_h(obj,pixel_side_h)
            if pixel_side_h<0
                error('pixel_side_h must be >0');
            end
            obj.pixel_side_h=pixel_side_h;
        end        
        
    end
end