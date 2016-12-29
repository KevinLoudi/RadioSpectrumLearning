function ab = alphabet(str)
% creates an alphabet object
%
% ab = alphabet(str)
%
% alphabet defines the "index" alphabet of a given string
% str - is a char-array 
%%%%
ab.ab_str = sort(unique(str));
ab = class(ab,'alphabet');
