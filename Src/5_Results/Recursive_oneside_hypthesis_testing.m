% ****************************************************
%  Note: recursive threshold to subtract signal part for noise/signal decision making
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

function [cut_point]=Recursive_oneside_hypthesis_testing(data, times,z_alph)
    if(nargin<2)
        error('Input parameters lacked!!!');
        exit;
    elseif(length(data)==0)
        error('Input data illegal!!!');
        exit;
    end
    tol_std=0.001;
    data_vec=data(:);
    %z_alph=1.645;%2.58;%1.645; % 95%confidence
    cut_point=zeros(1,times); %save threshold cut point
    d_mean=zeros(1,times); 
    d_std=zeros(1,times); 
    for ix=1:times
        d_mean(ix)=mean(data_vec);
        d_std(ix)=std(data_vec);
        cut_point(ix)=d_mean(ix)+z_alph*d_std(ix);
        remain_point=data_vec<cut_point(ix);
        %remain part
        data_vec=data_vec(remain_point);
        if ix>1 && abs(d_std(ix)-d_std(ix-1))<tol_std
                %display('reach the goal!!!');
                cut_point=cut_point(1:ix);
                break;
        end
    end
end