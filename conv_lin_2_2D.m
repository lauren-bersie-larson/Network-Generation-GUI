function [c] = conv_lin_2_2D(a, b)

% convert 1D array a ==> 2D array c with b cols 
%
% a -- row or column vector
% b -- scalar
% c -- 2D array with b cols
if size(a, 1) == 1 || size(a, 2) == 1
    c = reshape(a, b, [])';
else
    c = a;
end
