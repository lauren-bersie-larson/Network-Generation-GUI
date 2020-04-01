function [b] = conv_2D_2_lin(a)

% [b] = conv_2D_2_lin(a)
%
% convert 2D array a ==> 1D array b
%
% a -- 2D array
% b -- row vector
%
% last update -- sat aug 18 2012 -- mfh
if size(a,2) > 1
   b = reshape(a', 1, []);
else
    b = a;  
end
