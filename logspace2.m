function out = logspace2(a, b, n)
% This function is alternative version of logspace. 
% Whereas 'logspace' generates n log spaced points between decades 10^a 
% and 10^b, logspace2 generates n log spaced points between a and b.
%
% Inputs:
% a = interval start, scalar
% b = interval end, scalar
% n = number of points, integer
% 
% Ben Jancovich, 2022
% University of New South Wales
% This work is licensed under the Creative Commons Attribution 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
%
out = exp(linspace(log(a), log(b), n));

end