function translated_s = map(ab, s)
% maps string s to a corresponding string over alphanumeric-indices's alphabet 
% 
% translated_s = map(ab, s)
%
% s - is a sequence
% 
% Assumption:
% string "s" is alphanumeric.
%---
% Examples:
%>> ab = alphabet('abracadabra'); 
%>> s = map('abracadabra');
%>> sprintf('%d',s)
%
% ans =
%
% 12513141251
%
% Author: ron begleiter (http://www.cs.technion.ac.il/~ronbeg) 31 JULY 2007
%%%
for i=1:length(ab.ab_str)
    s(find(s==ab.ab_str(i))) = i; % NOTE: index of first char = "1" (in JAVA it should be "0") 
                                  %  thus, size(ab) = num of symbols+1 
end
translated_s = s;