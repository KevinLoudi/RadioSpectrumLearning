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

classdef stem_EM_options
    properties
        exit_toll=0.0001;                           %[double >0] (1x1) the EM algorithm stops if the relative norm between two consecutive iterations is below exit_toll
        max_iterations=100;                         %[integer >0](1x1) the EM algorithm stops if the number of iterations exceed max_iterations
        numeric_opt_type='single';                  %[string]    (1x1) 'single': then elements of the V_i matrices are numerically estimated one-by-one; 'full': the elements are jointly estimated.
        mstep_system_size=13500;                    %[integer >0](1x1) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
        compute_logL_at_all_steps=1;                %[boolean]   (1x1) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
        verbose=1;                                  %[boolean]   (1x1) 1: all the intermediate operations of the EM algorithm are displayed
        path_distributed_computing=[];              %[string]    (1x1) full or relative path of the folder to use for distributed computing
        timeout_distributed_computing=10000;        %[integer>0] (1x1) timeout in seconds when waiting for the data from the slaves
        timeout_node_search=10;                     %[integer>0] (1x1) timeout in seconds when looking for the available slaves
    end
    
    methods
        function obj = stem_EM_options(exit_toll,max_iterations,numeric_opt_type,mstep_system_size,compute_logL_at_all_steps,verbose,path_distributed_computing,timeout_distributed_computing,timeout_node_search)
            %DESCRIPTION: object constructor
            %
            %INPUT
            %<exit_toll>                     - [double >0]         (1x1) %(default: 0.0001) the EM algorithm stops if: 1) the relative norm of the model parameter vector between two iterations is below exit toll; 2) the relative norm of the observed data log-likelihood between two iterations is below exit toll (if computed)
            %<max_iterations>                - [integer >0]        (1x1) (default: 1000)  the EM algorithm stops if the number of iterations exceed max_iterations
            %<numeric_opt_type>              - [string]            (1x1) (default: 'single') 'single': then elements of the V_i matrices are numerically estimated one-by-one; 'full': the elements are jointly estimated.
            %<mstep_system_size>             - [integer >0]        (1x1) (default: 3500) if N_r(N_g)>mstep_system_size then theta_r and v_r (theta_g and v_g) are optimized by considering diagonal blocks of maximum dimension mstep_system_size
            %<compute_logL_at_all_steps>     - [boolean]           (1x1) (dafault: 0) 1: the observed data log-likelihood is evaluated at each iteration of the EM algorithm
            %<verbose>                       - [boolean]           (1x1) (default: 0) 1: all the intermediate operations of the EM algorithm are displayed
            %<path_distributed_computing>    - [string]            (1x1) (default: []) full or relative path of the folder to use for parallel computation
            %<timeout_distributed_computing> - [integer>0]         (1x1) (default: 10000) timeout in seconds when waiting for the data from the slaves 
            %<timeout_node_search>           - [integer>0]         (1x1) (default: 10) timeout in seconds when looking for the available slaves
            %
            %
            %OUTPUT
            %obj       - [stem_EM_options object] (1x1)
            
            if nargin>0
                obj.exit_toll=exit_toll;
            end
            if nargin>1
                obj.max_iterations=max_iterations;
            end
            if nargin>2
                obj.numeric_opt_type=numeric_opt_type;
            end
            if nargin>3
                obj.mstep_system_size=mstep_system_size;
            end
            if nargin>4
                obj.compute_logL_at_all_steps=compute_logL_at_all_steps;
            end
            if nargin>5
                obj.verbose=verbose;
            end
            if nargin>6
                obj.path_distributed_computing=path_distributed_computing;
            end
            if nargin>7
                obj.timeout_distributed_computing=timeout_distributed_computing;
            end
            if nargin>8
                obj.timeout_node_search=timeout_node_search;
            end
        end
        
        %Class set methods
        function obj = set.exit_toll(obj,exit_toll)
            if not(isempty(exit_toll))
                if exit_toll<=0
                    error('The exit_toll must be >0');
                end
                obj.exit_toll=exit_toll;
            end
        end
        
        function obj = set.max_iterations(obj,max_iterations)
            if not(isempty(max_iterations))
                if max_iterations<=0
                    error('max_iterations must be >0');
                end
                obj.max_iterations=max_iterations;
            end
        end
        
        function obj = set.numeric_opt_type(obj,numeric_opt_type)
            if not(isempty(numeric_opt_type))
                if not(strcmp(numeric_opt_type,'single'))&&not(strcmp(numeric_opt_type,'full'))
                    error('numeric_opt_type must be either ''single'' or ''full''');
                end
                obj.numeric_opt_type=numeric_opt_type;
            end
        end
        
        function obj = set.mstep_system_size(obj,mstep_system_size)
            if not(isempty(mstep_system_size))
                if mstep_system_size<=0
                    error('mstep_system_size must be >0');
                end
                obj.mstep_system_size=mstep_system_size;
            end
        end
        
        function obj = set.compute_logL_at_all_steps(obj,compute_logL_at_all_steps)
            if not(isempty(compute_logL_at_all_steps))
                if not(compute_logL_at_all_steps==0)&&not(compute_logL_at_all_steps==1)
                    error('compute_logL_at_all_steps must be either 0 or 1');
                end
                obj.compute_logL_at_all_steps=compute_logL_at_all_steps;
            end
        end
        
        function obj = set.verbose(obj,verbose)
            if not(isempty(verbose))
                if not(verbose==0)||not(verbose==1)
                    error('compute_logL_at_all_steps must be either 0 or 1');
                end
                obj.verbose=verbose;
            end
        end
        
        function obj = set.path_distributed_computing(obj,path_distributed_computing)
            if not(isempty(path_distributed_computing))
                if not(ischar(path_distributed_computing))
                    error('path_distributed_computing must be a string');
                else
                    obj.path_distributed_computing=path_distributed_computing;
                end
            end
        end
        
        function obj = set.timeout_distributed_computing(obj,timeout_distributed_computing)
            if not(isempty(timeout_distributed_computing))
                if timeout_distributed_computing<=0
                    error('timeout_distributed_computing must be >0');
                end
                obj.timeout_distributed_computing=timeout_distributed_computing;
            end
        end
        
        function obj = set.timeout_node_search(obj,timeout_node_search)
            if not(isempty(timeout_node_search))
                if timeout_node_search<=0
                    error('timeout_node_search must be >0');
                end
                obj.timeout_node_search=timeout_node_search;
            end
        end
          
    end
end