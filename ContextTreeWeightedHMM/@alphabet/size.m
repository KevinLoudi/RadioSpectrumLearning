function size_of_ab = size(ab)
% returns the size (=num of symbols+1) of the alphabet ab 
%
% size_of_ab = size(ab)
%
% NOTE: the size of the ab = num of symbols+1 due to MATLAB/JAVA
% implementation issues

size_of_ab = length(ab.ab_str)+1; % see NOTE