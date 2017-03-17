% Created on WED Mar 8th 13:25:45 2017
% Propose: Find grouped zeros in a seqence array
% Enviroment: Matlab 2015b
% @auththor: kevin

function [duration]=FindZerosBlock(sig)

 %threshold the vector to get a vector tsig of zeros
 tsig = (abs(sig) >= eps);
 %tsig = (abs(sig) > eps);
 
 %find the starting indices, ending indices, and duration of each string of zeroes
 dsig = diff([1 tsig 1]);
 startIndex = find(dsig < 0);
 endIndex = find(dsig > 0)-1;
 duration = endIndex-startIndex+1;
                     
end