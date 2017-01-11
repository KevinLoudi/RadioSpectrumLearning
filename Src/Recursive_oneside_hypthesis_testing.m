% ****************************************************
%  Note: recursive threshold to subtract signal part for noise/signal decision making
%  Author: Kevin
%  Date: 11th January, 2017
%  Environment: Matlab R2015b
%  Version: v1.0 (Last Modification Date: 11th January, 2017)
% ****************************************************

function Recursive_oneside_hypthesis_testing(data, times)
    [row,col]=size(data);
    data_vec=reshape(data,1,row*col);
    z_alph=1.645; % 95%confidence
    cut_point=zeros(1,times); %save threshold cut point
    for ix=1:times
        d_mean=mean(data_vec);
        d_std=std(data_vec);
        cut_point[ix]=d_mean+z_alph*d_std;
        %...............
    end
end